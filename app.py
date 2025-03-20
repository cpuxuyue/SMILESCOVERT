import streamlit as st
import subprocess
import sys

# 检查是否已安装rdkit
try:
    from rdkit import Chem
    from rdkit.Chem import Draw
except ImportError:
    st.error("""
    RDKit 未安装。请按照以下步骤安装：
    
    1. 在终端中运行：
    ```bash
    conda install -c conda-forge rdkit
    ```
    
    2. 安装完成后重新运行应用。
    """)
    st.stop()

from PIL import Image
import io

st.set_page_config(page_title="SMILES 结构查看器", layout="wide")

st.title("SMILES 结构查看器")
st.write("输入 SMILES 字符串，查看对应的化学结构")

# 创建两列布局
col1, col2 = st.columns(2)

with col1:
    smiles_input = st.text_input("输入 SMILES 字符串", "C1=CC=CC=C1")
    
    if smiles_input:
        try:
            # 创建 RDKit 分子对象
            mol = Chem.MolFromSmiles(smiles_input)
            if mol is None:
                st.error("无效的 SMILES 字符串")
            else:
                # 生成分子图片
                img = Draw.MolToImage(mol, size=(400, 400))
                
                # 将图片转换为字节流
                img_byte_arr = io.BytesIO()
                img.save(img_byte_arr, format='PNG')
                img_byte_arr = img_byte_arr.getvalue()
                
                # 显示图片
                st.image(img_byte_arr, caption="化学结构", use_container_width=True)
                
                # 添加下载按钮
                st.download_button(
                    label="下载结构图片",
                    data=img_byte_arr,
                    file_name="molecule.png",
                    mime="image/png"
                )
        except Exception as e:
            st.error(f"处理时出错: {str(e)}")

with col2:
    st.markdown("""
    ### 使用说明
    1. 在左侧输入框中输入 SMILES 字符串
    2. 系统会自动显示对应的化学结构
    3. 可以点击"下载结构图片"保存结构图
    
    ### 示例 SMILES
    - 苯: `C1=CC=CC=C1`
    - 乙醇: `CCO`
    - 阿司匹林: `CC(=O)OC1=CC=CC=C1C(=O)O`
    """) 