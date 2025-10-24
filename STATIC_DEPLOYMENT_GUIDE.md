# ClaimSpy UI - Static Deployment Guide

本指南帮助您将ClaimSpy UI从基于Python API服务器的动态版本转换为可部署到GitHub Pages的纯静态版本。

## 转换概述

### 原始架构 (动态版本)
- Python HTTP服务器 (`api_server.py`) 
- REST API端点用于数据访问
- 运行时文件读取和处理
- 需要Python环境

### 静态架构 (GitHub Pages版本)
- 纯静态HTML/CSS/JavaScript
- 预导出的JSON数据文件
- 基于配置的路径解析
- 无需服务器环境

## 使用步骤

### 1. 导出静态数据

运行导出脚本将动态数据转换为静态文件结构：

```bash
python3 export_static_data.py
```

这将：
- 扫描所有可用的评估文件夹
- 复制所有JSON数据文件到 `evaluation-interface/static-data/`
- 复制证据文件和模拟文件
- 生成配置文件和清单文件
- 创建根目录重定向页面

### 2. 本地测试

使用测试服务器在本地验证静态版本：

```bash
python3 test_static_server.py
```

然后访问：http://localhost:8000/evaluation-interface/

### 3. 配置模式切换

在 `evaluation-interface/static-config.js` 中：

```javascript
// 设置为true以启用静态模式（GitHub Pages）
isStaticMode: true

// 设置为false以使用API服务器模式（本地开发）
isStaticMode: false
```

### 4. 部署到GitHub Pages

1. 确保所有更改已提交到main分支
2. 推送到GitHub仓库
3. GitHub Actions将自动：
   - 运行静态数据导出
   - 构建Jekyll站点
   - 部署到GitHub Pages

## 文件结构

### 静态数据目录结构
```
evaluation-interface/static-data/
├── folders.json                          # 文件夹配置
├── alloys_0003/
│   ├── assessment.json                   # 评估数据
│   ├── trajectory.json                   # 轨迹数据
│   ├── mapping.json                      # 映射数据
│   ├── likert_score_prediction.json      # 预测数据
│   ├── evidence_source.json              # 证据源数据
│   ├── evidences/                        # 证据文件目录
│   │   ├── evidence1.pdf
│   │   └── evidence2.png
│   └── simulation-files/                 # 模拟文件目录
│       ├── manifest.json                 # 文件清单
│       ├── results.csv
│       ├── script.py
│       └── plot.png
└── computational_tools_0001/
    └── [similar structure]
```

### 配置文件
- `static-config.js` - 静态/动态模式配置
- `_config.yml` - Jekyll配置
- `.github/workflows/deploy.yml` - GitHub Actions工作流

## 主要更改

### JavaScript文件更改
所有API调用已更新为使用 `window.StaticConfig` 配置：

- `script.js` - 文件夹加载和数据加载
- `simulation_manager.js` - 模拟文件管理
- `pdf_controls.js` - PDF证据文件访问
- `evidence_manager.js` - 证据管理
- `code_explanation.js` - 代码解释文件访问
- `explanation_display.js` - 可视化图片显示

### 路径解析
- 动态模式：`/api/data/{folder}/{file}`
- 静态模式：`./static-data/{folder}/{file}`

## 优势

### 静态版本优势
✅ 可部署到GitHub Pages (免费托管)  
✅ 无需服务器维护  
✅ 快速加载和高可用性  
✅ 简单的部署流程  
✅ 版本控制和自动部署  

### 动态版本优势
✅ 运行时文件访问  
✅ 灵活的数据处理  
✅ 无需预处理步骤  
✅ 完整的文件系统访问  

## 切换模式

### 本地开发 (API服务器模式)
1. 设置 `isStaticMode: false` 在 `static-config.js`
2. 运行 `python3 api_server.py`
3. 访问 http://localhost:8080

### GitHub Pages部署 (静态模式)
1. 设置 `isStaticMode: true` 在 `static-config.js`
2. 运行 `python3 export_static_data.py`
3. 提交并推送到GitHub

## 疑难解答

### 常见问题

**问：文件找不到错误**  
答：确保运行了 `export_static_data.py` 并且文件存在于 `static-data/` 目录中

**问：模拟文件列表为空**  
答：检查 `manifest.json` 文件是否存在于相应的 `simulation-files/` 目录中

**问：证据文件无法加载**  
答：验证 `evidences/` 目录已正确复制到静态数据目录中

**问：GitHub Pages部署失败**  
答：检查GitHub Actions日志，确保所有必需的文件已提交到仓库

### 调试步骤
1. 检查浏览器控制台的错误信息
2. 验证网络请求的URL路径
3. 确认静态文件的完整性
4. 测试静态配置的开关状态