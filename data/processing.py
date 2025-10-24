import os
import json
from pathlib import Path
from openai import OpenAI

# ========== User-configurable section ==========
# Uncomment specific folders if you want to process only certain ones
folders = [
    # "alloys_0003",
    # "alloys_0004",
    # "alloys_0005",
    # Add more folders to process...
]


def auto_discover_folders(base_directory):
    """Automatically discover all folders that contain both assessment.json and trajectory.json."""
    discovered = []
    if base_directory.exists():
        for item in base_directory.iterdir():
            if item.is_dir():
                assessment_file = item / "assessment.json"
                trajectory_file = item / "trajectory.json"
                if assessment_file.exists() and trajectory_file.exists():
                    discovered.append(item.name)
    return sorted(discovered)


# Base paths and model setup
SCRIPT_DIR = Path(__file__).parent.absolute()
base_dir = SCRIPT_DIR / "sonnet45-full-run-codex-gpt5-1010" / "output" / "drop4"

# Process computational_tools_0001 for testing
folders = ["computational_tools_0001"]
model_name = "gpt-4o"
client = OpenAI(api_key=os.environ["OPENAI_TOKEN"])

# --- Prompt templates ---
ALIGN_PROMPT = """You are matching evidence from assessment to trajectory observations.
For each evidence item (with evidence ID), find which trajectory step(s) or reasoning section it corresponds to.

For simulation results: Match to step numbers where the simulation was performed (e.g., CALPHAD, CHGNet)
For reasoning: Match to the "Reasoning" section using -2, or to specific steps if the reasoning is step-specific

Return format: {"evidence_to_trajectory": {"ev_id1": [3, 4], "ev_id2": [-2], "ev_id3": [1]}} 
where:
- ev_id is the evidence identifier
- [3, 4] means steps 3 and 4
- [-2] means matches to the reasoning section
- [-1] means not found

Only match the evidence IDs provided in the input."""

DEPENDENCY_PROMPT = """You are analyzing step dependencies in a reasoning trajectory.
Each step has a thought that may reference or build upon previous observations.
Your task: for each step (0, 1, 2...), determine which previous step's observation it directly responds to or builds upon.

Important: Pay careful attention to:
1. Ordinal numbers in thoughts ("second search query" = step 1, "third search query" = step 2)
2. References to specific results within observations (when observations contain multiple results)
3. Content analysis to determine which specific result triggered the current thought

For dependencies, return:
- step_number: the step it depends on (-1 for step 0)
- result_index: if the observation contains multiple results (array), specify which result index (0-based) the step is responding to. Use null if not applicable or if observation is not an array.

Return JSON format: 
{
  "step_dependencies": {
    "0": {"depends_on_step": -1, "result_index": null},
    "1": {"depends_on_step": 0, "result_index": null},
    "2": {"depends_on_step": 1, "result_index": 0},
    ...
  }
}"""

# --- Main loop ---
for folder in folders:
    folder_path = base_dir / folder
    assessment_path = folder_path / "assessment.json"
    trajectory_path = folder_path / "trajectory.json"
    mapping_path = folder_path / "mapping.json"

    if not assessment_path.exists():
        print(f"  Missing assessment.json in {folder}")
        continue
    if not trajectory_path.exists():
        print(f"  Missing trajectory.json in {folder}")
        continue

    with open(assessment_path) as f:
        assessment = json.load(f)
    with open(trajectory_path) as f:
        trajectory = json.load(f)

    # --- Extract evidence and perform direct URL matching ---
    evidence_data = assessment.get("evidence", {})
    trajectory_data = trajectory.get("trajectory", {})

    # Classify evidence by type and handle direct observation mapping
    web_search_evidence = {}
    reasoning_simulation_evidence = {}
    evidence_to_trajectory = {}

    for ev_id, ev_info in evidence_data.items():
        ev_type = ev_info.get("type", "")
        source = ev_info.get("source", "")

        # Skip claim and other types - set to [-1]
        if ev_type in ["claim", "other"]:
            evidence_to_trajectory[ev_id] = [-1]
            continue

        # Handle direct observation mapping first
        if source.startswith("observation_") and source in trajectory_data:
            # Extract step number from observation key (e.g., "observation_0" -> 0)
            step_num = int(source.split("_")[1])
            evidence_to_trajectory[ev_id] = [step_num]
            print(
                f"    Direct observation mapping: {ev_id} -> {source} (step {step_num})")
            continue

        # Classify remaining evidence for further processing
        if ev_type == "web search":
            web_search_evidence[ev_id] = ev_info
        elif ev_type in ["reasoning", "simulation result"]:
            reasoning_simulation_evidence[ev_id] = ev_info

    # Direct URL matching for web search evidence
    for ev_id, ev_info in web_search_evidence.items():
        source_url = ev_info.get("source", "")
        if source_url:
            found = False
            for obs_key in trajectory_data.keys():
                if obs_key.startswith("observation_"):
                    obs_data = trajectory_data[obs_key]
                    if isinstance(obs_data, list):
                        for result_idx, result in enumerate(obs_data):
                            if isinstance(result, dict) and result.get("source") == source_url:
                                step_num = int(obs_key.split("_")[1])
                                evidence_to_trajectory[ev_id] = [
                                    step_num, result_idx]
                                found = True
                                break
                    if found:
                        break
            if not found:
                evidence_to_trajectory[ev_id] = [-1]  # Not found

    # Collect non-list observations for reasoning/simulation evidence
    non_list_observations = []
    for obs_key in sorted([k for k in trajectory_data.keys() if k.startswith("observation_")],
                          key=lambda x: int(x.split("_")[1])):
        obs_data = trajectory_data[obs_key]
        step_num = int(obs_key.split("_")[1])

        # Check if observation is not a list of dicts with source/citation
        if not isinstance(obs_data, list):
            non_list_observations.append(
                f"Step {step_num}: {json.dumps(obs_data, ensure_ascii=False)[:300]}...")
        else:
            # Check if it's not a list of search results
            if not all(isinstance(item, dict) and "source" in item and "citation" in item for item in obs_data):
                non_list_observations.append(
                    f"Step {step_num}: {json.dumps(obs_data, ensure_ascii=False)[:300]}...")

    # Add reasoning section if present
    if "reasoning" in trajectory_data:
        reasoning_text = trajectory_data["reasoning"]
        non_list_observations.append(f"Reasoning: {reasoning_text[:500]}...")

    # Prepare evidence for LLM matching (only remaining reasoning and simulation results)
    evidence_for_llm = []
    for ev_id, ev_info in reasoning_simulation_evidence.items():
        evidence_text = f"Evidence ID: {ev_id}\nType: {ev_info.get('type')}\nSource: {ev_info.get('source')}\nCitation: {ev_info.get('citation', '')[:200]}..."
        evidence_for_llm.append(evidence_text)

    print(f"  Processing {folder} ...")
    print(f"    Found {len(web_search_evidence)} web search evidence items")
    print(
        f"    Found {len(reasoning_simulation_evidence)} reasoning/simulation evidence items")
    print(
        f"    Direct observation matches: {len([k for k, v in evidence_to_trajectory.items() if len(v) == 1 and v[0] >= 0])}")
    print(
        f"    Direct URL matches: {len([k for k, v in evidence_to_trajectory.items() if len(v) == 2])}")
    print(
        f"    Evidence set to [-1]: {len([k for k, v in evidence_to_trajectory.items() if v == [-1]])}")

    # Only call LLM if there are evidence items that need LLM processing
    if evidence_for_llm:
        # --- Construct model input for evidence matching ---
        user_input = f"""
Evidence to match (reasoning and simulation results only):
{chr(10).join(evidence_for_llm)}

Non-list trajectory observations and reasoning:
{chr(10).join(non_list_observations)}

{ALIGN_PROMPT}
"""

        # --- Call structured-output model for explanation to trajectory mapping ---
        completion1 = client.chat.completions.create(
            model=model_name,
            messages=[
                {"role": "system",
                    "content": "You are a precise reasoning alignment assistant."},
                {"role": "user", "content": user_input}
            ],
            response_format={
                "type": "json_schema",
                "json_schema": {
                    "name": "EvidenceToTrajectory",
                    "schema": {
                        "type": "object",
                        "properties": {
                            "evidence_to_trajectory": {
                                "type": "object",
                                "additionalProperties": {
                                    "type": "array",
                                    "items": {"type": "integer"}
                                }
                            }
                        },
                        "required": ["evidence_to_trajectory"]
                    }
                }
            },
            temperature=0
        )

        explanation_mapping = completion1.choices[0].message.content
        print(f"    LLM evidence mapping: {explanation_mapping}")

        # Parse and merge LLM results
        evidence_json = json.loads(explanation_mapping)
        if "evidence_to_trajectory" in evidence_json:
            evidence_to_trajectory.update(
                evidence_json["evidence_to_trajectory"])
    else:
        print(f"    No evidence items require LLM processing")
        evidence_json = {"evidence_to_trajectory": {}}

    # --- Extract step thoughts and observations for dependency analysis ---
    step_thoughts = []
    step_keys = sorted(
        [k for k in trajectory_data.keys() if k.startswith("thought_")],
        key=lambda x: int(x.split("_")[1])
    )

    for key in step_keys:
        step_idx = int(key.split("_")[1])
        thought_content = trajectory_data[key]

        # Include previous observations for context
        prev_observations = []
        for i in range(step_idx):
            obs_key = f"observation_{i}"
            if obs_key in trajectory_data:
                obs_data = trajectory_data[obs_key]
                if isinstance(obs_data, list):
                    obs_summary = f"Observation {i} (array with {len(obs_data)} results)"
                    for j, result in enumerate(obs_data):
                        if isinstance(result, dict) and 'source' in result:
                            obs_summary += f"\n  [{j}] {result.get('source', 'Unknown source')}: {result.get('citation', 'No citation')[:100]}..."
                        else:
                            obs_summary += f"\n  [{j}] {str(result)[:100]}..."
                else:
                    obs_summary = f"Observation {i}: {str(obs_data)[:200]}..."
                prev_observations.append(obs_summary)

        context = "\n".join(
            prev_observations) if prev_observations else "No previous observations"
        step_thoughts.append(
            f"Step {step_idx}:\nThought: {thought_content}\nPrevious context:\n{context}\n---")

    dependency_input = f"""
Step thoughts from trajectory:
{chr(10).join(step_thoughts)}

{DEPENDENCY_PROMPT}
"""

    print(f"    Analyzing step dependencies...")

    # --- Call structured-output model for step dependency analysis ---
    completion2 = client.chat.completions.create(
        model=model_name,
        messages=[
            {"role": "system",
                "content": "You are a precise step dependency analysis assistant."},
            {"role": "user", "content": dependency_input}
        ],
        response_format={
            "type": "json_schema",
            "json_schema": {
                "name": "StepDependencies",
                "schema": {
                    "type": "object",
                    "properties": {
                        "step_dependencies": {
                            "type": "object",
                            "additionalProperties": {
                                "type": "object",
                                "properties": {
                                    "depends_on_step": {"type": "integer"},
                                    "result_index": {"type": ["integer", "null"]}
                                },
                                "required": ["depends_on_step", "result_index"]
                            }
                        }
                    },
                    "required": ["step_dependencies"]
                }
            }
        },
        temperature=0
    )

    dependency_mapping = completion2.choices[0].message.content
    print(f"    Step dependency mapping: {dependency_mapping}")

    # --- Combine both mappings into final result ---
    dependency_json = json.loads(dependency_mapping)

    # Evidence mapping is already complete in evidence_to_trajectory
    combined_mapping = {
        "evidence_to_trajectory": evidence_to_trajectory,
        "step_dependencies": dependency_json["step_dependencies"]
    }

    # --- Save result ---
    with open(mapping_path, "w") as f:
        json.dump(combined_mapping, f, indent=2)

    print(f"  Saved combined mapping.json to {folder}")
