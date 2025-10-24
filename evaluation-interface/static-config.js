// Static Data Configuration for GitHub Pages deployment
// This file manages the transition from API-based dynamic loading to static file loading

window.StaticConfig = {
    // Set to true for GitHub Pages deployment, false for local development with API server
    isStaticMode: true,

    // Base path for static data files
    staticDataPath: './static-data',

    // API endpoints (used when isStaticMode is false)
    apiEndpoints: {
        folders: '/api/folders',
        data: '/api/data',
        simulationFiles: '/api/simulation-files',
        simulationFile: '/api/simulation-file'
    },

    // Static file paths (used when isStaticMode is true)
    staticPaths: {
        folders: './static-data/folders.json',
        data: (folder, filename) => `./static-data/${folder}/${filename}`,
        evidences: (folder, filename) => `./static-data/${folder}/evidences/${filename}`,
        simulationFiles: (folder) => `./static-data/${folder}/simulation-files`,
        simulationFile: (folder, filename) => `./static-data/${folder}/simulation-files/${filename}`
    },

    // Utility functions for getting correct paths
    getFoldersUrl() {
        return this.isStaticMode ? this.staticPaths.folders : this.apiEndpoints.folders;
    },

    getDataUrl(folder, filename) {
        return this.isStaticMode ?
            this.staticPaths.data(folder, filename) :
            `${this.apiEndpoints.data}/${folder}/${filename}`;
    },

    getEvidenceUrl(folder, filename) {
        return this.isStaticMode ?
            this.staticPaths.evidences(folder, filename) :
            `${this.apiEndpoints.data}/${folder}/evidences/${filename}`;
    },

    getSimulationFilesUrl(runType, folder) {
        // In static mode, we don't need the runType as files are organized by folder
        return this.isStaticMode ?
            this.staticPaths.simulationFiles(folder) :
            `${this.apiEndpoints.simulationFiles}/${runType}/${folder}`;
    },

    getSimulationFileUrl(runType, folder, filename) {
        // In static mode, we don't need the runType as files are organized by folder
        return this.isStaticMode ?
            this.staticPaths.simulationFile(folder, filename) :
            `${this.apiEndpoints.simulationFile}/${runType}/${folder}/${filename}`;
    },

    // Check if a file exists (for static mode)
    async fileExists(url) {
        if (!this.isStaticMode) {
            return true; // API server will handle 404s
        }

        try {
            const response = await fetch(url, { method: 'HEAD' });
            return response.ok;
        } catch (error) {
            return false;
        }
    },

    // Load simulation files list for static mode
    async loadSimulationFilesList(folder) {
        if (!this.isStaticMode) {
            throw new Error('loadSimulationFilesList should only be used in static mode');
        }

        try {
            // Load the manifest file created during export
            const manifestUrl = `./static-data/${folder}/simulation-files/manifest.json`;
            const response = await fetch(manifestUrl);

            if (response.ok) {
                const manifest = await response.json();
                return { files: manifest.files || [] };
            } else {
                console.warn(`No simulation files manifest found for ${folder}`);
                return { files: [] };
            }
        } catch (error) {
            console.error(`Error loading simulation files manifest for ${folder}:`, error);
            return { files: [] };
        }
    }
};