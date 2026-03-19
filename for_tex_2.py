import re
from sympy import symbols, simplify, solve, Eq

def extract_vars(expr_str):
    """从字符串中提取所有连续字母（变量名）"""
    return set(re.findall(r'[a-zA-Z]+', expr_str))

def parse_user_substitutions():
    print("请输入变量替换规则（格式：var = expression），每行一条，输入空行结束：")
    print("示例：a = c + 2")
    print("     b = d + 3\n")
    
    rules = {}
    all_vars = set()
    while True:
        line = input(">>> ").strip()
        if not line:
            break
        if '=' not in line:
            print("⚠️  格式错误，应为 'var = expression'，跳过。")
            continue
        var, expr = line.split('=', 1)
        var = var.strip()
        expr = expr.strip()
        if not var.isalpha():
            print(f"⚠️  变量名 '{var}' 应为纯字母，跳过。")
            continue
        rules[var] = expr
        all_vars.add(var)
        all_vars.update(extract_vars(expr))
    
    # 基础变量 = 出现过但不在等号左边的变量
    lhs_vars = set(rules.keys())
    base_vars = sorted(all_vars - lhs_vars)  # 排序保证顺序一致
    
    print(f"\n检测到基础变量（自由参数）: {base_vars}")
    return rules, base_vars


def build_normalizer_advanced(rules, base_vars):
    # 创建符号
    base_symbols = {name: symbols(name) for name in base_vars}
    all_symbols = base_symbols.copy()
    
    # 解析替换规则为sympy表达式
    rule_exprs = {}
    for var, expr_str in rules.items():
        try:
            # 创建这个变量的符号
            var_sym = symbols(var)
            all_symbols[var] = var_sym
            
            # 解析表达式
            expr = eval(expr_str, {"__builtins__": {}}, base_symbols)
            rule_exprs[var_sym] = simplify(expr)
        except Exception:
            print(f"警告：无法解析规则 {var} = {expr_str}")
    
    def normalize_index(idx_str):
        idx_str = idx_str.strip()
        if idx_str.isdigit():
            return idx_str
        
        try:
            # 解析输入表达式
            expr = eval(idx_str, {"__builtins__": {}}, all_symbols)
            simplified = simplify(expr)
            
            # 应用所有已知的替换规则
            for var_sym, var_expr in rule_exprs.items():
                simplified = simplified.subs(var_sym, var_expr)
            
            # 简化
            simplified = simplify(simplified)
            
            # 转换为字符串
            result = str(simplified).replace(" ", "")
            return result
            
        except Exception as e:
            # 如果无法解析，返回原始字符串（无空格）
            return idx_str.replace(" ", "")
    
    return normalize_index



def build_normalizer(rules, base_vars):
    # 动态创建符号
    sym_dict = {name: symbols(name) for name in base_vars}
    
    # 构建命名空间
    local_ns = sym_dict.copy()
    
    # 解析规则（多次迭代处理依赖）
    parsed = {}
    remaining = rules.copy()
    max_iter = 10
    for _ in range(max_iter):
        progress = False
        to_remove = []
        for var, expr_str in remaining.items():
            try:
                expr = eval(expr_str, {"__builtins__": {}}, local_ns)
                parsed[var] = simplify(expr)
                local_ns[var] = parsed[var]
                to_remove.append(var)
                progress = True
            except Exception:
                continue
        for var in to_remove:
            del remaining[var]
        if not remaining:
            break
        if not progress:
            break
    
    if remaining:
        print(f"⚠️  以下规则无法解析（依赖缺失或语法错误）：{list(remaining.keys())}")

    def normalize_index(idx_str):
        idx_str = idx_str.strip()
        if idx_str.isdigit():
            return idx_str
        if idx_str in parsed:
            return str(simplify(parsed[idx_str])).replace(" ", "")
        # 尝试作为表达式直接解析（使用基础符号）
        try:
            expr = eval(idx_str, {"__builtins__": {}}, sym_dict)
            return str(simplify(expr)).replace(" ", "")
        except Exception:
            return idx_str.replace(" ", "")  # 保留原样

    return normalize_index

# --- 以下函数与之前相同，仅需传入 normalize_func ---

def normalize_m_label(m_str, normalize_func):
    if not m_str.startswith("M_{"):
        return m_str
    match = re.match(r"M_\{(.*)\}", m_str)
    if not match:
        return m_str
    inner = match.group(1)
    if '|' in inner:
        left, right = inner.split('|', 1)
        left_parts = [p.strip() for p in left.split(',')]
        right_parts = [p.strip() for p in right.split(',')]
        norm_left = [normalize_func(p) for p in left_parts]
        norm_right = [normalize_func(p) for p in right_parts]
        normalized_inner = ",".join(norm_left) + "|" + ",".join(norm_right)
    else:
        parts = [p.strip() for p in inner.split(',')]
        normalized_inner = ",".join(normalize_func(p) for p in parts)
    return f"M_{{{normalized_inner}}}"

def read_blocks(filepath):
    with open(filepath, 'r') as f:
        content = f.read()
    blocks = [block.strip() for block in content.split('\n\n') if block.strip()]
    return blocks

def parse_file_normalized(filepath, normalize_func):
    blocks = read_blocks(filepath)
    if len(blocks) < 5:
        raise ValueError(f"文件 {filepath} 至少需要5个块（#1 到 #5）")
    u_lines = [line.strip() for line in blocks[4].split('\n') if line.strip()]
    n = len(u_lines)
    if len(blocks) < 5 + n:
        raise ValueError(f"文件 {filepath} 缺少数据块（需要 {n} 个，实际 {len(blocks)-5}）")
    data_blocks = blocks[5:5+n]
    
    u_norm = [normalize_m_label(m, normalize_func) for m in u_lines]
    data_map = {}
    for i, m_orig in enumerate(u_lines):
        key = normalize_m_label(m_orig, normalize_func)
        data_lines = [line.strip() for line in data_blocks[i].split('\n') if line.strip()]
        data_norm = [normalize_m_label(d, normalize_func) for d in data_lines]
        data_map[key] = data_norm
    return u_norm, data_map

def is_contained(A_file, B_file, normalize_func):
    uA, dataA = parse_file_normalized(A_file, normalize_func)
    uB, dataB = parse_file_normalized(B_file, normalize_func)
    
    setA, setB = set(uA), set(uB)
#    if not setA.issubset(setB):
    if not setA == setB:
        print(f"问题在开始处这里:")
        print(f"setA-setB:{setA-setB}")
        print(f"setB-setA:{setB-setA}")
        return False
    for m in setA:
        if m not in dataB:
            return False
        if not set(dataA[m]).issubset(set(dataB[m])):
            print(f"问题在这里{m}:")
            print(f"{set(dataA[m])-set(dataB[m])}")
            print("")
            return False
    return True

def for_tex_contained(input1,input2):
    # 1. 获取用户规则和基础变量
    rules, base_vars = parse_user_substitutions()
    
    # 2. 构建标准化函数
    normalize_func = build_normalizer_advanced(rules, base_vars)
    print(f"***:{rules}")
    
    # 3. 输入文件
#    fileA = input("\n请输入文件 A 的路径: ").strip()
#    fileB = input("请输入文件 B 的路径: ").strip()
    fileA = input1
    fileB = input2
    # 4. 对比
    try:
        result = is_contained("test//"+fileA+".txt", "test//"+fileB+".txt", normalize_func)
        print(f"\n✅ 结果: B {'包含' if result else '不包含'} A")
    except Exception as e:
        print(f"\n❌ 错误: {e}")

# ===== 主程序 =====
#if __name__ == "__main__":
#    for_tex_contained("input2","input3")
