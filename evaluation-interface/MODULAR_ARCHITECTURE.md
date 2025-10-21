# ClaimSpy UI - Modular Architecture

## Module Structure

The ClaimSpy UI has been refactored into a modular architecture for better maintainability and code organization. Here's the breakdown of the modules:

### Core Modules

#### 1. `script.js` - Main Application Module
- **Purpose**: Main application logic, folder management, data loading, and UI coordination
- **Key Functions**:
  - `loadAvailableFolders()` - Load available assessment folders
  - `loadSelectedFolder()` - Load and display selected assessment data
  - `loadAssessmentData()` - Parse and display assessment results
  - `displayClaim()`, `displayClaimCard()`, `displaySystemScores()` - UI display functions
  - Panel management (JSON panel, trajectory flow coordination)

#### 2. `evidence_manager.js` - Evidence Management Module
- **Purpose**: Handle all evidence-related functionality including simulation files and local evidence
- **Key Functions**:
  - `loadSimulationFiles()` - Load and manage simulation files
  - `showSimulationFile()` - Display simulation files with visualizations
  - `viewLocalEvidence()` - Handle local PDF evidence viewing
  - `formatCsvAsTable()` - Format CSV files as HTML tables
  - `getEvidenceTypeClass()` - Styling for different evidence types
- **Exports**: `window.evidenceManager` object with all public functions

#### 3. `explanation_display.js` - Explanation Display Module  
- **Purpose**: Handle the display and interaction with explanations and their associated evidence
- **Key Functions**:
  - `displayExplanations()` - Main function to render all explanations
  - `formatTrajectoryText()` - Format trajectory text with markdown-like styling
  - `formatSearchResults()` - Format search results as structured cards
  - `toggleExplanationContent()`, `markExplanation()` - Explanation interaction
- **Exports**: `window.explanationDisplay` object with all public functions

#### 4. `simulation_manager.js` - Simulation Management Module
- **Purpose**: Handle all simulation file loading, display, and processing
- **Key Functions**:
  - `loadSimulationFiles()` - Load simulation files from API
  - `showSimulationFile()` - Display simulation files with visualizations
  - `formatCsvAsTable()` - Format CSV files as HTML tables
  - `formatFileSize()` - Format file sizes for display
  - `parseCsvLine()` - Parse CSV content
- **Exports**: `window.simulationManager` object with all simulation functions

#### 5. `code_explanation.js` - Code Explanation Module (Pre-existing)
- **Purpose**: Handle Python code files with AI-generated explanations
- **Key Functions**:  
  - `isExplanationFile()` - Identify explanation files
  - `showSimulationFileWithExplanation()` - Display Python files with explanations
  - `loadExplanationData()` - Load explanation data for Python files
- **Exports**: `window.codeExplanation` object + backward compatibility globals

### Supporting Modules

#### 5. `trajectory_flow.js` - Trajectory Visualization (Pre-existing)
- **Purpose**: Handle trajectory flow visualization in left panel
- **Exports**: `window.trajectoryFlow` object

#### 6. `pdf_controls.js` - PDF Viewer (Pre-existing)  
- **Purpose**: Enhanced PDF viewing with PDF.js integration
- **Key Functions**: PDF rendering, highlighting, navigation

#### 7. `panel_resize.js` - Panel Management (Pre-existing)
- **Purpose**: Handle resizable panels and UI interactions

## Module Loading Order

The modules are loaded in this specific order in `index.html`:

```html
### Module Loading Order

```html
<!-- Essential order for proper initialization -->
<script src="evidence_manager.js"></script>    <!-- Evidence management -->
<script src="explanation_display.js"></script> <!-- Explanation display -->  
<script src="simulation_manager.js"></script>  <!-- Simulation file handling -->
<script src="code_explanation.js"></script>    <!-- Code explanation handling -->
<script src="script.js"></script>              <!-- Main application -->
```
```

## Inter-Module Communication

### Dependencies
- `explanation_display.js` depends on `simulationManager` for simulation file operations and `evidenceManager` for evidence data
- `evidence_manager.js` delegates simulation-related functions to `simulationManager` 
- `evidence_manager.js` optionally uses `code_explanation.js` for identifying explanation files
- `script.js` coordinates with all modules through their exported objects

### Global Objects
Each module exports its public API through a global `window` object:
- `window.evidenceManager` - Evidence management functions
- `window.explanationDisplay` - Explanation display functions
- `window.simulationManager` - Simulation file handling functions
- `window.codeExplanation` - Code explanation functions
- `window.trajectoryFlow` - Trajectory visualization

### Module API Usage Examples

```javascript
// Evidence module usage
window.evidenceManager.viewLocalEvidence(filepath, mimetype);
window.evidenceManager.getEvidenceTypeClass(filepath);

// Explanation module usage  
window.explanationDisplay.displayExplanations(data);
window.explanationDisplay.formatTrajectoryText(text);

// Simulation module usage
window.simulationManager.loadSimulationFiles(problemFolder, runType);
window.simulationManager.showSimulationFile(data, filename);
window.simulationManager.formatCsvAsTable(csvContent);

// Code explanation module usage
window.codeExplanation.isExplanationFile(filename);
window.codeExplanation.showSimulationFileWithExplanation(content, filename);
```

### Wrapper Functions
`script.js` provides wrapper functions for backward compatibility:
- `showSimulationFile()` - Delegates to `simulationManager.showSimulationFile()`
- `toggleEvidenceSection()` - Delegates to `evidenceManager.toggleEvidenceSection()`
- `viewLocalEvidence()` - Delegates to `evidenceManager.viewLocalEvidence()`
- `toggleExplanationContent()` - Delegates to `explanationDisplay.toggleExplanationContent()`

## Benefits of Modular Architecture

1. **Separation of Concerns**: Each module has a specific responsibility
2. **Maintainability**: Easier to locate and modify specific functionality
3. **Reusability**: Modules can be reused across different parts of the application
4. **Testing**: Individual modules can be tested in isolation
5. **Code Organization**: Related functions are grouped together logically
6. **Reduced File Size**: Main script.js is significantly smaller and focused

## Migration Notes

The refactoring maintains full backward compatibility. All existing functionality continues to work as before, but the underlying implementation is now properly modularized.

Key changes:
- Evidence management functions moved from `script.js` to `evidence_manager.js`
- Explanation display logic moved from `script.js` to `explanation_display.js`
- Simulation file handling moved from `evidence_manager.js` to `simulation_manager.js`
- Global variables for evidence data now managed through `evidenceManager.setEvidenceData()`
- Function calls now use proper module delegation where appropriate
- Simulation-related functions now properly separated from evidence management