# Scientific Claim Assessment Interface

## Overview
A comprehensive web-based interface for evaluating scientific claims with AI-powered analysis, trajectory visualization, and evidence extraction. The interface supports modular architecture, real-time content analysis, and integrated Chrome extension functionality for enhanced research workflows.

## Key Features

### 🔧 **Modular Architecture**
- **Separated codebase**: Clean separation of HTML, CSS, and JavaScript files
- **`index.html`**: Streamlined HTML structure with external resource links
- **`styles.css`**: Complete CSS styling for responsive three-panel layout
- **`script.js`**: Full JavaScript functionality including AI integration
- **Maintainable code**: Easy to edit, debug, and collaborate on

### 📊 **Data Management & Visualization**
- **Folder-based loading**: Dynamic selection from organized result directories
- **Multiple data formats**: Support for assessment, trajectory, mapping, and prediction files
- **JSON context panel**: Interactive navigation and data highlighting
- **Real-time updates**: Automatic loading and refresh capabilities

### 🤖 **AI-Powered Analysis**
- **Content scraping**: Multi-proxy web content extraction with fallback strategies
- **OpenAI integration**: GPT-4 powered content analysis and passage identification
- **Relevance scoring**: Automated assessment of content relevance to claims
- **Access restriction handling**: Graceful fallback for paywalled/restricted content

### 🌐 **Chrome Extension Integration**
- **Direct page highlighting**: Real-time text highlighting in source documents
- **Cross-tab communication**: Seamless integration between interface and browser tabs
- **Enhanced analysis**: Full-page content analysis with visual feedback

### 📈 **Trajectory Visualization**
- **Step-by-step reasoning**: Visual representation of AI reasoning processes
- **Color-coded elements**: Thoughts (blue), Tools (green), Observations (orange)
- **Interactive exploration**: Expandable sections with detailed information
- **Mapping integration**: Links between explanations and reasoning steps

### 📚 **Evidence Management**
- **Multi-source evidence**: Support for web searches, academic papers, and databases
- **Citation formatting**: Proper academic citation display
- **Content caching**: Intelligent caching to avoid redundant web requests
- **Reading mode**: Focused content display with relevant passage highlighting

## Data Structure
The interface expects data to be organized in the following structure:
```
data/sonnet45-full-run-codex-gpt5-1010/output/dry-run/
├── alloys_0001/
│   ├── assessment.json      # Main assessment results
│   ├── trajectory.json      # AI reasoning steps
│   ├── mapping.json         # Links explanations to trajectory
│   └── likert_score_prediction.json  # Scoring predictions
├── alloys_0002/
│   ├── assessment.json
│   ├── trajectory.json
│   ├── mapping.json
│   └── likert_score_prediction.json
└── ...
```

## Quick Start

### 1. Setup Configuration
```bash
# Copy the template to create your local configuration
cd evaluation-interface
cp local-config.template.js local-config.js

# Edit local-config.js to add your OpenAI API key (optional)
# window.OPENAI_TOKEN = 'your-openai-api-key-here';
```

### 2. Start the API Server
```bash
# From the evaluation-interface directory
python3 api_server.py
```

### 3. Access the Interface
- Open http://localhost:8080 in your browser
- The modular interface will load with separated CSS/JS files

### 4. Configure API Keys
- Click "Settings" to configure OpenAI API key for content analysis
- Keys are stored locally and remembered across sessions

### 5. Select and Analyze Results
- Use the dropdown to select a result folder (e.g., alloys_0003)
- View assessment scores, explanations, and evidence
- Explore trajectories with color-coded reasoning steps
- Analyze evidence sources with AI-powered content extraction

### 6. Chrome Extension (Optional)
- Install the ClaimSpy Chrome extension for enhanced functionality
- Enable direct page highlighting and cross-tab analysis
- Get real-time content analysis in source document tabs

## API Reference

### Core Endpoints

#### GET /api/folders
Returns list of available result folders
```json
{
  "folders": ["alloys_0001", "alloys_0002", "alloys_0003", ...]
}
```

#### GET /api/data/{folder}/{filename}
Returns JSON content of specified file
- **folder**: Result folder name (e.g., alloys_0003)
- **filename**: One of: assessment.json, trajectory.json, mapping.json, likert_score_prediction.json

#### Static File Serving
- **GET /**: Serves index.html
- **GET /styles.css**: CSS styling
- **GET /script.js**: JavaScript functionality
- **GET /config.js**: Configuration settings
- **GET /local-config.js**: Local configuration overrides

## Data File Formats

### assessment.json
```json
{
  "type": "assessment",
  "problem_id": "alloys_0003",
  "likert_score": -1,
  "continuous_score": 0.08,
  "confidence": 0.62,
  "explanation": [
    {
      "text": "Explanation text...",
      "evidence": ["ev_id1", "ev_id2"]
    }
  ],
  "evidence": {
    "ev_id1": {
      "type": "web search",
      "source": "https://...",
      "citation": "Citation text..."
    }
  }
}
```

### trajectory.json
```json
{
  "trajectory": {
    "thought_1": "Analysis text...",
    "tool_name_1": "search",
    "tool_call_1": { "query": "..." },
    "observation_1": "Search results...",
    "thought_2": "Further analysis...",
    "tool_name_2": "calphad",
    "tool_call_2": { "composition": "..." },
    "observation_2": "Simulation results..."
  }
}
```

### mapping.json
```json
{
  "explanation_to_trajectory": {
    "1": {
      "matched_thoughts": [1],
      "matched_tools": ["search"],
      "matched_observations": [1]
    },
    "2": {
      "matched_thoughts": [2, 3],
      "matched_tools": ["calphad"],
      "matched_observations": [2, 3]
    }
  }
}
```

## Features

### Visual Elements
- **Folder selection dropdown**: Easy navigation between results
- **Color-coded trajectory types**: Thoughts (blue), Tools (green), Observations (orange)
- **Collapsible sections**: Expandable trajectory and evidence details
- **Formatted displays**: Clean presentation of all data types

### Integration
- **Seamless loading**: All files loaded automatically when folder is selected
- **Evidence analysis**: Full AI-powered content analysis with OpenAI integration
- **JSON context panel**: Navigate and highlight specific data elements
- **Responsive design**: Works on desktop and mobile devices

## Configuration

### Server Configuration
Update paths in `api_server.py`:
```python
DATA_BASE_PATH = '/path/to/your/data/directory'
STATIC_BASE_PATH = '/path/to/evaluation-interface'
PORT = 8080
```

### OpenAI API Setup
1. Get API key from OpenAI
2. Configure in interface settings
3. Key is stored locally in browser

### Problem Claims
Add custom claims in `script.js` `getProblemClaimForFolder()` function:
```javascript
const problemClaims = {
    'new_domain_001': 'Your custom claim text here',
    // Add more as needed
};
```

### Chrome Extension
1. Load extension from `/Claimspy-reader-extension/` directory
2. Enable developer mode in Chrome
3. Configure extension permissions

## Advanced Features

### Content Analysis Strategies
- **Multi-proxy scraping**: Automatic fallback through multiple CORS proxies
- **Content extraction**: Smart text extraction from various webpage formats
- **Access restriction handling**: Graceful degradation for restricted content
- **Caching system**: Intelligent caching to avoid redundant requests

### Trajectory Analysis
- **Step mapping**: Links explanations to specific reasoning steps
- **Color coding**: Visual distinction between thoughts, tools, and observations
- **Interactive exploration**: Expandable sections with detailed information

### Evidence Integration
- **Multi-source support**: Web searches, academic papers, databases
- **Citation management**: Proper academic citation formatting
- **Relevance scoring**: AI-powered assessment of content relevance

## Troubleshooting

### Common Issues
- **Modular loading errors**: Ensure all CSS/JS files are in same directory
- **No folders appear**: Check DATA_BASE_PATH in api_server.py
- **AI analysis fails**: Verify OpenAI API key configuration
- **Chrome extension issues**: Check extension permissions and popup blockers
- **Content scraping fails**: Try different proxy services or direct access

### Development Mode
- **Demo data**: Interface falls back to simulated data if API server unavailable
- **Console logging**: Detailed error messages in browser console
- **Cache clearing**: Use browser dev tools to clear evidence cache

## File Architecture

### Core Files
- **`index.html`**: Main interface structure
- **`styles.css`**: Complete styling (740+ lines)
- **`script.js`**: Full functionality (1500+ lines)
- **`api_server.py`**: Backend API server

### Supporting Files
- **`config.js`**: Base configuration
- **`local-config.template.js`**: Template for local configuration
- **`local-config.js`**: Local overrides (not tracked in git)
- **`README.md`**: This documentation
- **`.gitignore`**: Prevents committing sensitive files

### Chrome Extension
- **`Claimspy-reader-extension/`**: Complete browser extension
- **Cross-tab communication**: Seamless integration
- **Page highlighting**: Real-time text highlighting

## Development

### Adding New Features
1. **CSS changes**: Edit `styles.css` for styling
2. **JavaScript functionality**: Add to `script.js`
3. **HTML structure**: Modify `index.html` sparingly
4. **Server endpoints**: Extend `api_server.py` as needed

### Code Organization
- **Modular structure**: Clean separation of concerns
- **Maintainable codebase**: Easy to debug and extend
- **Professional architecture**: Industry-standard file organization

## Configuration & Security

### Path Configuration
The application now uses relative paths automatically:
- **Data path**: `../data/sonnet45-full-run-codex-gpt5-1010/output/dry-run/` (relative to script)
- **Static path**: Current directory of `api_server.py`
- **No hardcoded user paths**: Safe for version control and distribution

### API Key Management
```bash
# 1. Copy template to create local config
cp local-config.template.js local-config.js

# 2. Edit local-config.js to add your API key
# window.OPENAI_TOKEN = 'sk-proj-your-key-here';

# 3. local-config.js is automatically ignored by git
```

### Security Best Practices
- **API keys**: Never commit to version control
- **Local config**: Use `local-config.js` for sensitive settings
- **Path validation**: Server validates all file paths
- **CORS enabled**: Secure cross-origin handling
- **File restrictions**: Only allowed data files are served