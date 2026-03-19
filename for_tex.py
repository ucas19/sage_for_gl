def read_blocks(filepath):
    with open(filepath, 'r') as f:
        content = f.read()
    # 按空行分割，过滤空块
    blocks = [block.strip() for block in content.split('\n\n') if block.strip()]
    return blocks

def format_expression_lines(lines_str):
    """格式化 #3, #4：第一行无 +，其余每行前加 + 并换行"""
    if not lines_str.strip():
        return ""
    lines = [line.strip() for line in lines_str.strip().split('\n') if line.strip()]
    if not lines:
        return ""
    result = lines[0]
    for line in lines[1:]:
        result += "\n+ " + line
    return result

def format_u_set(u_lines_str):
    """格式化 #5：每行末尾加逗号，最后一行不加"""
    lines = [line.strip() for line in u_lines_str.strip().split('\n') if line.strip()]
    if not lines:
        return ""
    # 所有行加逗号，然后去掉最后一个逗号
    with_commas = [line + "," for line in lines]
    with_commas[-1] = lines[-1]  # 最后一行不加逗号
    return "\n".join(with_commas)

def m_to_p_in_dollar(m_expr):
    """将 M_{...} 转为 $P_{...}$"""
    if m_expr.startswith("M_{"):
        p_expr = "P_" + m_expr[2:]  # replace M_ with P_
        return f"${p_expr}$"
    else:
        # 如果不是 M_ 开头，原样加美元符号？
        return f"${m_expr}$"

def format_proj_covers(u_lines_str):
    """将 U 集合的每一行转为 $P_{...}$，然后每4个一组，用 $a,b,c,d$ \\ 形式"""
    lines = [line.strip() for line in u_lines_str.strip().split('\n') if line.strip()]
    if not lines:
        return []
    
    # 转为 $P_{...}$
    p_list = [m_to_p_in_dollar(line) for line in lines]
    
    # 每4个一组
    return p_list

def format_cell_content(cell_str):
    """将 cell 内容中的每行 M_{...} 转为 $M_{...}$，并保持原分行"""
    lines = [line.strip() for line in cell_str.strip().split('\n') if line.strip()]
    if not lines:
        return ""
    # 每行加 $...$
    dollar_lines = [f"${line}$" for line in lines]
    # 每4行加 \\ 换行？但样例是手动分行的，我们保持原样，只加 $

    grouped = []
    for i in range(0, len(dollar_lines), 4):
        group = dollar_lines[i:i+4]
        group[0] = " " + group[0]
        grouped.append(",\n ".join(group))
#    return ",\n".join(dollar_lines)
    return grouped

def generate_tex_from_file(filepath):
    blocks = read_blocks(filepath)
    
    # 至少需要 5 个 blocks
    if len(blocks) < 5:
        raise ValueError("Input file must contain at least 5 blocks (separated by blank lines).")
    
    block0 = blocks[1]  # #1
    block1 = blocks[0]  # #2
    block2 = blocks[2]  # #3
    block3 = blocks[3]  # #4
    block4 = blocks[4]  # #5
    
    # 获取 U 行数
    u_lines = [line.strip() for line in block4.strip().split('\n') if line.strip()]
    n = len(u_lines)
    
    # 检查是否有足够的 cell blocks
    if len(blocks) < 5 + n:
        raise ValueError(f"Expected {n} cell blocks after block4, but got {len(blocks) - 5}.")
    
    cell_blocks = blocks[5:5+n]  # #7, #9, ...
    
    # 格式化各部分
#    p_mu = block0.strip()
    v_index_2 = block1.strip()  # 用于 V_{#2}
    v_index_4 = block3.strip()  # 用于 V_{#4}（虽然应相同）
    v_index_4 = ','.join(v_index_4) # 用于 V_{#4}（虽然应相同）
    
    p_mu = format_expression_lines(block0)
    expr_line1 = format_expression_lines(block2)
    expr_line2 = format_expression_lines(block3)  # 注意：#4 是简化式，可能含 \sum
    
    u_formatted = format_u_set(block4)
    
    # 生成 Proj. Cover 列表（#6, #8, ...）
    proj_covers = format_proj_covers(block4)  # 每个元素是一串 "$a,b,c,d$ \\"
    
    # 生成 cell 内容（#7, #9, ...）
    cell_contents = [format_cell_content(cb) for cb in cell_blocks]
    
    # 构建表格行
    table_rows = []
    for i in range(n):
        proj = proj_covers[i] if i < len(proj_covers) else "$ $"
        cell = cell_contents[i]
        row = f"""$ {proj.replace(' \\\\', '')} $ & \\makecell{{\n{cell}\n}}\\\\ \n"""
        # 但注意：proj_covers[i] 已经是 "$a,b,c,d$ \\"，我们只需取内容
        # 更简单：proj_covers 是列表 of string like "$P1, P2$", 不带 \\
        # 重新调整 format_proj_covers
        pass
    
    # 重新定义：proj_items 是 [$P1$, $P2$, ...]
    u_lines_clean = [line.strip() for line in block4.strip().split('\n') if line.strip()]
    proj_items = [m_to_p_in_dollar(line) for line in u_lines_clean]
    
    # 现在构建表格行
    table_body = ""
    for i in range(n):
        proj = proj_items[i]
        cell = cell_contents[i]
        cell_add = ""
        for cell_tem in cell:
            cell_add +=cell_tem
            cell_add += ", \\\\" + "\n" 
        cell_add = cell_add[:-5] + "\n"

        table_body += f""" {proj} & \\makecell{{\n{cell_add}}}\\\\ \n"""
        table_body += " & \\\\ \n"
    table_body = table_body[:-8]
    
    # LaTeX 模板（使用 .format，注意转义 {{}}）
    tex_template = r"""        Let $\mu:=\lambda + (0,0,0,0|0) $. 
        \begin{{equation}}\notag
          \begin{{aligned}}
            P_\mu = {p_mu}
          \end{{aligned}}
        \end{{equation}}
        \begin{{equation}}\notag
          \begin{{aligned}}
            Pr_{{\lambda}}(P_\mu \otimes V_{{{v_index_2}}}) &= 
{expr_line1} \\
                                               &=
{expr_line2}
          \end{{aligned}}
        \end{{equation}}
      Similarly, the set of indeterminate composition factors is $U = \{{
{u_formatted}
\}}$. We have $Pr_{{\lambda}}(P_\mu \otimes V_{{{v_index_2}}}) = P_\lambda  $ by \cref{{tab:verma-analysis_41_a_0}}. 
\begin{{table}}[htbp]
\centering
\caption{{Calculation of Verma modules in $\operatorname{{Filt}}(P_*) \setminus \operatorname{{Filt}}(\prl(P_\lambda \otimes V_{{{v_index_2}}})) $ }} \label{{tab:verma-analysis_41_a_0}}
\begin{{tabular}}{{cc}}
\toprule
\textbf{{Proj. Cover}} & \textbf{{$\operatorname{{Filt}}(P_*) \setminus \operatorname{{Filt}}(\prl(P_\mu \otimes V_{{{v_index_2}}})) $}}  \\
     & \textbf{{(via Prop.～1)}}  \\
\midrule
{table_body}
\bottomrule
\end{{tabular}}
\vspace{{0.5em}}
\parbox{{\linewidth}}{{\small
\textbf{{Interpretation:}} 
This table displays only selected computational results. This partial presentation does not affect our conclusion.
}}
\end{{table}}"""
    
    return tex_template.format(
        p_mu=p_mu,
        v_index_2=v_index_2,
        v_index_4=v_index_4,
        expr_line1=expr_line1,
        expr_line2=expr_line2,
        u_formatted=u_formatted,
        table_body=table_body
    )

# 主程序
def for_tex_entry(adress):
    tex_output = generate_tex_from_file("test//"+adress+".txt")
    print(tex_output)


#--------------------------------------------------------------

