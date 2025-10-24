# ClaimSpy UI - Scientific Claim Assessment Interface

## 🎯 Overview

ClaimSpy UI is a comprehensive web-based interface for evaluating scientific claims with AI-powered analysis, trajectory visualization, and evidence extraction. The system provides interactive tools for researchers to analyze scientific assessments, explore reasoning trajectories, and examine supporting evidence.

## 🚀 Deployment Options

ClaimSpy UI now supports **two deployment modes**:

### 📱 **Static Mode (GitHub Pages)** - *Recommended for sharing results*
- ✅ **Zero setup** - Works directly on GitHub Pages
- ✅ **No server required** - Pure HTML/CSS/JavaScript  
- ✅ **Fast loading** - All data pre-exported as JSON files
- ✅ **Easy sharing** - Public URL for research presentations

### 🔧 **Dynamic Mode (Local Server)** - *Best for development*
- ✅ **Live data access** - Real-time file system integration
- ✅ **Development flexibility** - Easy data updates and modifications
- ✅ **Full functionality** - Complete API-based data access

### Quick Start

**For GitHub Pages deployment (Static Mode):**
```bash
# 1. Export static data
python3 export_static_data.py

# 2. Test locally (optional)
python3 test_static_server.py
# Visit: http://localhost:8000/evaluation-interface/

# 3. Push to GitHub - automatic deployment via GitHub Actions
```

**For local development (Dynamic Mode):**
```bash
# Set isStaticMode: false in static-config.js, then:
python3 api_server.py
# Visit: http://localhost:8080
```

📖 **Detailed Guide**: See [STATIC_DEPLOYMENT_GUIDE.md](STATIC_DEPLOYMENT_GUIDE.md) for complete instructions.



## ✨ Key Features## Changes Made



### 📊 **Interactive Assessment Dashboard**### 1. API Server (api_server.py)

- **Claim Evaluation**: View detailed assessments of scientific claims with Likert scores, continuous scores, and confidence metrics**Before:**

- **Dynamic Folder Selection**: Browse and switch between different assessment results```python

- **Real-time Loading**: Automatic detection and loading of available assessment foldersDATA_BASE_PATH = '/Users/litaoyan/Documents/Research/HAILMAIER-C/ClaimSpy_UI/data/sonnet45-full-run-codex-gpt5-1010/output/dry-run'

STATIC_BASE_PATH = '/Users/litaoyan/Documents/Research/HAILMAIER-C/ClaimSpy_UI/evaluation-interface'

### 🔍 **Evidence Analysis & Visualization**```

- **Evidence Explorer**: Interactive examination of supporting evidence with source links and citations

- **AI-Powered Content Analysis**: Automatic extraction and highlighting of relevant passages from scientific documents**After:**

- **Multi-source Evidence**: Support for web searches, academic papers, simulations, and experimental data```python

# Use relative paths based on the script location

### 🧠 **Trajectory Visualization**SCRIPT_DIR = Path(__file__).parent.absolute()

- **Reasoning Path Display**: Step-by-step visualization of AI reasoning processesPROJECT_ROOT = SCRIPT_DIR.parent

- **Interactive Trajectory**: Expandable sections showing thoughts, tool usage, and observationsDATA_BASE_PATH = PROJECT_ROOT / 'data' / 'sonnet45-full-run-codex-gpt5-1010' / 'output' / 'dry-run'

- **Color-coded Elements**: Visual distinction between thoughts (blue), tools (green), and observations (orange)STATIC_BASE_PATH = SCRIPT_DIR

- **Mapping Integration**: Links between explanations and specific reasoning steps```



### 🤖 **AI Integration**### 2. Frontend JavaScript (script.js)

- **OpenAI GPT Integration**: Optional AI-powered content analysis and passage identification**Before:**

- **Smart Caching**: Intelligent caching system to avoid redundant API calls```javascript

- **Fallback Strategies**: Graceful handling of restricted or inaccessible contentconst DATA_BASE_PATH = '/Users/litaoyan/Documents/Research/HAILMAIER-C/final outputs/sonnet45-full-run-codex-gpt5-1010/output/dry-run';

```

### 📱 **Modern Web Interface**

- **Three-panel Layout**: Efficient use of screen space with collapsible panels**After:**

- **Responsive Design**: Works across different screen sizes and devices```javascript

- **JSON Context Panel**: Interactive JSON viewer with syntax highlighting// Data directory configuration

- **Reading Mode**: Focused content display with relevant passage highlighting// Note: This path is now handled by the API server, so this constant is kept for reference only

```

## 🚀 Quick Start

### 3. API Key Security (local-config.js)

### 1. Prerequisites**Before:**

- Python 3.7+ ```javascript

- Modern web browserwindow.OPENAI_TOKEN = 'sk-proj-Krafg7PtrY8J3Em8vIQNT3BlbkFJhhReko6hYIRRHPtL5bdg';

- (Optional) OpenAI API key for enhanced AI features```



### 2. Installation**After:**

```bash```javascript

# Clone the repository// window.OPENAI_TOKEN = 'your-openai-api-key-here';

git clone https://github.com/YanLitao/ClaimSpy_UI.git```

cd ClaimSpy_UI

### 4. New Files Created

# Navigate to the interface directory

cd evaluation-interface#### .gitignore

```- Prevents committing sensitive files (local-config.js, API keys, etc.)

- Excludes common temporary and system files

### 3. Setup Configuration (Optional)

```bash#### local-config.template.js

# Copy the configuration template- Template file for users to copy and configure their own API keys

cp local-config.template.js local-config.js- Safe to commit to version control



# Edit local-config.js to add your OpenAI API key (optional)## Benefits

# window.OPENAI_TOKEN = 'your-openai-api-key-here';

```### Security

- ✅ No hardcoded personal paths

### 4. Start the Server- ✅ No API keys in version control

```bash- ✅ Template-based configuration system

# From the evaluation-interface directory- ✅ Proper .gitignore protection

python3 api_server.py

```### Portability

- ✅ Works on any system without modification

### 5. Access the Interface- ✅ Automatic path detection based on script location

- Open your browser and navigate to: **http://localhost:8080**- ✅ Relative path structure maintained

- The interface will automatically load available assessment folders

### Maintainability

## 📁 Data Structure- ✅ Clear separation of public and private configuration

- ✅ Easy setup process for new users

The system expects data to be organized in the following structure:- ✅ Professional project structure



```## Setup for New Users

ClaimSpy_UI/

├── evaluation-interface/           # Web interface1. **Clone the repository**

│   ├── api_server.py              # Backend API server   ```bash

│   ├── index.html                 # Main interface   git clone <repository-url>

│   ├── script.js                  # Frontend logic   cd ClaimSpy_UI

│   ├── styles.css                 # Styling   ```

│   └── local-config.js            # User configuration

├── data/                          # Data directory2. **Setup local configuration**

│   └── sonnet45-full-run-codex-gpt5-1010/   ```bash

│       └── output/   cd evaluation-interface

│           └── dry-run/           # Assessment results   cp local-config.template.js local-config.js

│               ├── alloys_0001/   # Example assessment folder   # Edit local-config.js to add API key if needed

│               │   ├── assessment.json      # ⭐ Main assessment results   ```

│               │   ├── trajectory.json      # 🧠 AI reasoning steps  

│               │   ├── mapping.json         # 🔗 Links explanations to trajectory3. **Run the server**

│               │   └── likert_score_prediction.json  # 📊 Scoring predictions   ```bash

│               ├── alloys_0002/   python3 api_server.py

│               ├── batteries_0001/   ```

│               ├── superconductors_0001/

│               └── ...4. **Access the interface**

└── README.md   - Open http://localhost:8080

```

## Path Structure

### 🔑 Key Files in Each Assessment Folder```

ClaimSpy_UI/

#### `assessment.json` ⭐├── evaluation-interface/          # Frontend and API server

Contains the main assessment results:│   ├── api_server.py             # Main server (uses relative paths)

```json│   ├── index.html

{│   ├── script.js

  "type": "assessment",│   ├── styles.css

  "problem_id": "alloys_0001",│   ├── config.js

  "likert_score": -1,│   ├── local-config.template.js   # Template (safe to commit)

  "continuous_score": 0.08,│   ├── local-config.js           # User config (ignored by git)

  "confidence": 0.62,│   └── README.md

  "explanation": [├── data/                         # Data directory (relative path)

    {│   └── sonnet45-full-run-codex-gpt5-1010/

      "text": "Detailed explanation of the assessment...",│       └── output/

      "evidence": ["ev_id1", "ev_id2"]│           └── dry-run/

    }│               ├── alloys_0001/

  ],│               ├── alloys_0002/

  "evidence": {│               └── ...

    "ev_id1": {└── .gitignore                    # Protects sensitive files

      "type": "web search",```

      "source": "https://...",

      "citation": "Citation text..."## Verification

    }- [x] No hardcoded user paths remain

  }- [x] API keys removed from tracked files

}- [x] Relative paths work correctly

```- [x] .gitignore protects sensitive data

- [x] Template system for configuration

#### `trajectory.json` 🧠- [x] README updated with new instructions

Contains AI reasoning steps:- [x] Code passes syntax validation

```json

{The project is now ready for GitHub upload and safe for public sharing!
  "trajectory": {
    "thought_1": "Initial analysis...",
    "tool_name_1": "search",
    "tool_call_1": {"query": "..."},
    "observation_1": "Search results...",
    "thought_2": "Further analysis...",
    "observation_2": "Additional findings..."
  }
}
```

#### `mapping.json` 🔗
Links explanations to trajectory steps:
```json
{
  "explanation_to_trajectory": {
    "1": 1,  // Explanation 1 → Trajectory step 1
    "2": 3,  // Explanation 2 → Trajectory step 3
    "3": 4   // Explanation 3 → Trajectory step 4
  }
}
```

## 🎮 How to Use

### Basic Navigation
1. **Select Assessment**: Use the dropdown menu to choose an assessment folder
2. **View Results**: Examine the claim, scores, and explanations in the main panel
3. **Explore Evidence**: Click on evidence items to view source details
4. **Check Trajectory**: Use "Show Trajectory" buttons to see reasoning steps

### Advanced Features
- **AI Analysis**: Click the document icon (📄) next to evidence to analyze source content
- **JSON Context**: Use the left panel to explore raw JSON data structure
- **Reading Mode**: Click evidence sources to open in reading mode with highlighted passages
- **Cross-referencing**: Navigate between explanations and their corresponding trajectory steps

### Panel Controls
- **Left Panel**: JSON context viewer (toggle with arrow button)
- **Center Panel**: Main assessment content
- **Right Panel**: Evidence and reading mode (toggle with arrow button)

## 🔧 Configuration Options

### API Server Configuration
- **Port**: Default 8080 (modify in `api_server.py`)
- **Data Path**: Automatically detected relative to project structure
- **CORS**: Enabled for development

### AI Features Configuration
- **OpenAI Integration**: Optional, requires API key in `local-config.js`
- **Content Analysis**: Automatic passage extraction and relevance scoring
- **Caching**: Smart caching to minimize API usage

## 🛠️ Development & Data Processing

### Processing New Data
The project includes a data processing script for generating mappings:

```bash
cd data
export OPENAI_API_KEY="your-key-here"
python3 processing.py
```

This script:
- Processes assessment and trajectory files
- Generates mapping.json files using AI analysis
- Supports batch processing of multiple folders

### Adding New Assessments
1. Create a new folder in `data/sonnet45-full-run-codex-gpt5-1010/output/dry-run/`
2. Add required JSON files (`assessment.json`, `trajectory.json`)
3. Optionally run `processing.py` to generate `mapping.json`
4. Refresh the web interface to see the new assessment

## 🚨 Troubleshooting

### Common Issues
- **No folders showing**: Check that data directory exists and contains valid assessment folders
- **File not found errors**: Ensure each folder has required JSON files
- **AI features not working**: Verify OpenAI API key is set in `local-config.js`
- **Server won't start**: Check that port 8080 is available

### API Endpoints
- `GET /api/folders` - List available assessment folders
- `GET /api/data/{folder}/{file}` - Retrieve specific data files  
- `GET /api/health` - Server health check

## 🤝 Contributing
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Test with sample data
5. Submit a pull request

## 📄 License
[Add your license information here]

---

**Ready to start analyzing scientific claims?** Launch the server and explore the interactive assessment interface! 🚀