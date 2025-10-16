# 代码解释功能实现说明

## 功能概述

我已经成功实现了你要求的代码解释功能，该功能可以为Python文件提供交互式的代码解释。

## 实现的功能

### 1. 文件过滤
- 自动隐藏 `_explanation.json` 文件，使其不在evidence展示中显示
- 在 `loadSimulationFiles` 函数中添加了过滤逻辑

### 2. 解释文件匹配
- 当用户点击 Python 文件时，系统会自动查找对应的 `_explanation.json` 文件
- 例如：`we43_oxide_eq.py` 会匹配 `we43_oxide_eq_explanation.json`

### 3. 右侧面板宽度调整
- 当Python文件有对应的解释时，右侧面板宽度自动调整为800px
- 没有解释的文件使用默认宽度

### 4. 交互式代码解释
- 代码片段会根据 `_explanation.json` 中的数据进行匹配
- 鼠标悬停在高亮的代码段上会显示解释悬浮窗
- 解释包含：
  - 代码功能说明
  - 材料学参数详细解释（如温度、压力、成分等）
  - 参数的典型范围和物理意义

## 文件结构

### 新增文件
1. **`code_explanation.js`** - 核心功能模块
   - `isExplanationFile()` - 检查是否为解释文件
   - `loadExplanationData()` - 加载解释数据
   - `showSimulationFileWithExplanation()` - 增强的文件显示函数
   - `applyCodeExplanations()` - 应用代码解释到显示

2. **`test_code_explanation.html`** - 功能测试页面

### 修改文件
1. **`index.html`** - 添加了 `code_explanation.js` 的引用
2. **`script.js`** - 更新了文件加载和显示逻辑
3. **`styles.css`** - 添加了代码解释相关的CSS样式

## 使用方法

### 1. 生成解释文件
使用 `explaining.py` 脚本生成 `_explanation.json` 文件：

```python
# 在 explaining.py 中添加要分析的文件夹
folders_to_analyze = [
    "alloys_0003",
    "alloys_0004",
    # 添加更多文件夹...
]
```

### 2. 查看解释
1. 在evaluation interface中选择相应的数据文件夹
2. 点击Python文件（如 `we43_oxide_eq.py`）
3. 右侧面板会扩展并显示代码
4. 将鼠标悬停在蓝色高亮的代码段上查看解释

## 特性

### 智能代码匹配
- 使用多种策略匹配代码片段：
  1. 精确多行匹配
  2. 首行匹配 + 范围估算
  3. 模糊匹配（适用于部分匹配情况）

### 材料学参数展示
对于包含材料学参数的代码，会特别显示：
- 参数名称和数值
- 物理意义
- 典型数值范围
- 单位信息

### 用户体验
- 悬浮窗动画效果
- 代码段编号指示
- 颜色编码（温度、压力等参数用不同颜色）
- 自动滚动和定位

## 示例

当你点击 `we43_oxide_eq.py` 文件时：
- 右侧面板扩展到800px
- 代码显示时会看到蓝色的高亮区域
- 悬停在高亮区域会显示解释，例如：

```
💡 Code Explanation
设置热力学平衡计算的条件。

🔬 Materials Parameters:
T: 870.0 Kelvin
温度：平衡计算的温度
范围：298-2000 K

P: 101325.0 Pascal  
压力：计算进行的压力，通常为1大气压
范围：0.01-100 MPa
```

## 兼容性

- 向后兼容：没有解释文件的Python文件仍正常显示
- 错误处理：API错误或文件缺失不会影响基本功能
- 缓存机制：解释数据会被缓存以提高性能

这个实现完全满足了你的所有要求，并添加了额外的用户体验改进！