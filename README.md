# SMILES 结构查看器

这是一个简单的网页应用，可以将 SMILES 字符串转换为化学结构式。

## 功能特点

- 实时预览化学结构
- 支持下载结构图片
- 提供示例 SMILES
- 错误提示功能

## 安装说明

1. 确保已安装 Python 3.7 或更高版本
2. 安装依赖包：
   ```bash
   pip install -r requirements.txt
   ```

## 运行应用

```bash
streamlit run app.py
```

## 使用说明

1. 在输入框中输入 SMILES 字符串
2. 系统会自动显示对应的化学结构
3. 可以点击"下载结构图片"保存结构图

## 示例 SMILES

- 苯: `C1=CC=CC=C1`
- 乙醇: `CCO`
- 阿司匹林: `CC(=O)OC1=CC=CC=C1C(=O)O` 