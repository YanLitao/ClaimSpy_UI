import os
import json
from pathlib import Path
from openai import OpenAI

# ========== User-configurable section ==========
folders = [
    # "alloys_0003",
    #"alloys_0004",
    "alloys_0005",
    # "batteries_0001",
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

# Uncomment this line if you prefer auto-discovery instead of manual folder list
# folders = auto_discover_folders(base_dir)


# Base paths and model setup
SCRIPT_DIR = Path(__file__).parent.absolute()
base_dir = SCRIPT_DIR / "sonnet45-full-run-codex-gpt5-1010" / "output" / "dry-run"
model_name = "gpt-4o"
client = OpenAI(api_key=os.environ["OPENAI_TOKEN"])

# --- Prompt template ---
ALIGN_PROMPT = """You are given two JSON files: assessment.json and trajectory.json.
Each assessment contains multiple explanations, and each trajectory contains multiple thought/observation pairs.
Your task: for each explanation (indexed from 1), find the single most relevant observation number (0-based, as shown)
that best supports the explanation content.
Return only a JSON object mapping explanation index (1-based) → observation index (0-based)."""

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

    # --- Extract explanations and observations ---
    explanations = [
        f"{i+1}. {exp['text']}" for i, exp in enumerate(assessment.get("explanation", []))
    ]

    observations = []
    trajectory_data = trajectory.get("trajectory", {})
    observation_keys = sorted(
        [k for k in trajectory_data.keys() if k.startswith("observation_")],
        key=lambda x: int(x.split("_")[1])
    )
    for key in observation_keys:
        idx = int(key.split("_")[1])
        val = trajectory_data[key]
        obs_text = val if isinstance(
            val, str) else json.dumps(val, ensure_ascii=False)
        observations.append(f"{idx}. {obs_text}")

    # --- Construct model input ---
    user_input = f"""
Assessment explanations:
{chr(10).join(explanations)}

Trajectory observations:
{chr(10).join(observations)}

{ALIGN_PROMPT}
"""

    print(f"  Processing {folder} ...")
    print(f"    Found {len(explanations)} explanations")
    print(f"    Found {len(observations)} observations")

    # --- Call structured-output model ---
    completion = client.chat.completions.create(
        model=model_name,
        messages=[
            {"role": "system", "content": "You are a precise reasoning alignment assistant."},
            {"role": "user", "content": user_input}
        ],
        response_format={
            "type": "json_schema",
            "json_schema": {
                "name": "ExplanationToTrajectory",
                "schema": {
                    "type": "object",
                    "properties": {
                        "explanation_to_trajectory": {
                            "type": "object",
                            "additionalProperties": {"type": "integer"}
                        }
                    },
                    "required": ["explanation_to_trajectory"]
                }
            }
        },
        temperature=0
    )

    # --- Extract structured output using beta API ---
    mapping = completion.choices[0].message.content
    print(f"    Model output: {mapping}")

    # --- Save result ---
    with open(mapping_path, "w") as f:
        f.write(mapping)

    print(f"  Saved structured mapping.json to {folder}")
