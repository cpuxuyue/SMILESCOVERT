import streamlit as st
import subprocess
import sys
import pandas as pd
import os
import base64

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
import zipfile

st.set_page_config(page_title="SMILES 结构查看器", layout="wide")

st.title("SMILES 结构查看器")
st.write("输入 SMILES 字符串或上传 CSV 文件，查看对应的化学结构")

# 创建两个标签页
tab1, tab2 = st.tabs(["单个转换", "批量转换"])

with tab1:
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

with tab2:
    st.markdown("""
    ### 批量转换说明
    1. 上传包含 SMILES 列的 CSV 文件
    2. 选择包含 SMILES 的列名
    3. 系统会自动生成并显示所有结构的图片
    4. 可以下载单个结构图片
    """)
    
    uploaded_file = st.file_uploader("上传 CSV 文件", type=['csv'])
    
    if uploaded_file is not None:
        try:
            # 读取CSV文件
            df = pd.read_csv(uploaded_file)
            
            # 显示列名选择器
            smiles_column = st.selectbox("选择包含 SMILES 的列", df.columns)
            
            if smiles_column:
                # 创建两列布局用于显示图片
                col1, col2 = st.columns(2)
                
                # 处理每个SMILES
                success_count = 0
                error_count = 0
                error_smiles = []
                
                for idx, row in df.iterrows():
                    smiles = str(row[smiles_column])
                    try:
                        mol = Chem.MolFromSmiles(smiles)
                        if mol is not None:
                            # 生成分子图片
                            img = Draw.MolToImage(mol, size=(400, 400))
                            
                            # 将图片转换为字节流
                            img_byte_arr = io.BytesIO()
                            img.save(img_byte_arr, format='PNG')
                            img_byte_arr = img_byte_arr.getvalue()
                            
                            # 在对应的列中显示图片
                            if idx % 2 == 0:
                                with col1:
                                    st.image(img_byte_arr, caption=f"SMILES: {smiles}", use_container_width=True)
                                    st.download_button(
                                        label=f"下载结构 {idx+1}",
                                        data=img_byte_arr,
                                        file_name=f"molecule_{idx+1}.png",
                                        mime="image/png",
                                        key=f"download_{idx}"
                                    )
                            else:
                                with col2:
                                    st.image(img_byte_arr, caption=f"SMILES: {smiles}", use_container_width=True)
                                    st.download_button(
                                        label=f"下载结构 {idx+1}",
                                        data=img_byte_arr,
                                        file_name=f"molecule_{idx+1}.png",
                                        mime="image/png",
                                        key=f"download_{idx}"
                                    )
                            
                            success_count += 1
                        else:
                            error_count += 1
                            error_smiles.append(smiles)
                    except Exception as e:
                        error_count += 1
                        error_smiles.append(smiles)
                
                # 显示统计信息
                st.success(f"处理完成！成功: {success_count}, 失败: {error_count}")
                
                # 如果有失败的SMILES，显示错误信息
                if error_smiles:
                    with st.expander("查看失败的 SMILES"):
                        st.write("以下 SMILES 处理失败：")
                        for smiles in error_smiles:
                            st.code(smiles)
                
        except Exception as e:
            st.error(f"处理CSV文件时出错: {str(e)}") 