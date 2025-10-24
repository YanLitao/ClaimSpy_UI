import os
import json
from pathlib import Path
from openai import OpenAI


def explain_python_code():
    """
    Analyze Python files in specified folders and generate explanations using OpenAI GPT-4o.
    """

    # 1. List of folder names to analyze
    folders_to_analyze = [
        # "sandbox-dry-run/alloys_0003",
        # Add more folder names here as needed
        "sandbox-drop4/computational_tools_0001",
    ]

    # 2. Base directory path
    SCRIPT_DIR = Path(__file__).parent.absolute()
    base_dir = SCRIPT_DIR / "sonnet45-full-run-codex-gpt5-1010"

    # Initialize OpenAI client
    client = OpenAI(api_key=os.environ["OPENAI_TOKEN"])

    # 4. Prompt template for GPT-4o
    prompt_template = """
You are a materials science and computational simulation expert.  
Your task is to **split the entire Python code below into several meaningful snippets**.  
The snippets must:
1. **Cover every line of code in order, with no omissions.**  
   Each line must appear exactly once in one snippet.  
2. **Keep the original code unaltered** (preserve spacing and indentation).  
3. For each snippet, explain its purpose clearly.  
4. If it includes materials science parameters (temperature, pressure, composition, potentials, etc.), explain their physical meaning and typical range.

Output must strictly follow the JSON schema below.

Code:
{code_content}
"""

    # 5. Structured output schema
    response_format = {
        "type": "json_schema",
        "json_schema": {
            "name": "materials_code_analysis",
            "schema": {
                "type": "object",
                "properties": {
                    "snippets": {
                        "type": "array",
                        "items": {
                            "type": "object",
                            "properties": {
                                "code": {"type": "string"},
                                "explanation": {"type": "string"},
                                "materials_parameters": {
                                    "type": "array",
                                    "items": {
                                        "type": "object",
                                        "properties": {
                                            "name": {"type": "string"},
                                            "value": {"type": "string"},
                                            "units": {"type": "string"},
                                            "meaning": {"type": "string"},
                                            "typical_range": {"type": "string"}
                                        },
                                        "required": ["name", "value", "units", "meaning", "typical_range"],
                                        "additionalProperties": False
                                    },
                                    "default": []
                                }
                            },
                            "required": ["code", "explanation", "materials_parameters"],
                            "additionalProperties": False
                        }
                    }
                },
                "required": ["snippets"],
                "additionalProperties": False
            },
            "strict": True
        }
    }

    # Process each folder
    for folder_name in folders_to_analyze:
        folder_path = base_dir / folder_name
        print(f"Processing folder: {folder_name}")

        if not folder_path.exists():
            print(f"Warning: Folder {folder_path} does not exist. Skipping...")
            continue

        # 3. Find all .py files in the folder
        py_files = list(folder_path.glob("*.py"))

        if not py_files:
            print(f"No Python files found in {folder_path}")
            continue

        for py_file in py_files:
            print(f"  Analyzing file: {py_file.name}")

            try:
                # Read the Python file
                with open(py_file, 'r', encoding='utf-8') as f:
                    code_content = f.read()

                # Prepare the prompt
                prompt = prompt_template.format(code_content=code_content)

                # Call OpenAI API
                response = client.chat.completions.create(
                    model="gpt-4o",
                    messages=[
                        {
                            "role": "user",
                            "content": prompt
                        }
                    ],
                    response_format=response_format,
                    temperature=0  # Lower temperature for more consistent output
                )

                # Parse the response
                analysis_result = json.loads(
                    response.choices[0].message.content)

                # 6. Save the output
                output_filename = f"{py_file.stem}_explanation.json"
                output_path = folder_path / output_filename

                with open(output_path, 'w', encoding='utf-8') as f:
                    json.dump(analysis_result, f, indent=4, ensure_ascii=False)

                print(f"    Saved explanation to: {output_filename}")

            except Exception as e:
                print(f"    Error processing {py_file.name}: {str(e)}")
                continue

    print("Analysis complete!")


if __name__ == "__main__":
    explain_python_code()
