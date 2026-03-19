import ast
import re

def parse_simple_config(filename):
    """
    解析简单的配置文件格式
    格式示例：
    # 注释
    setA = [7]
    setB = 1,2,3,4,5,6,7,8,9,10,11,12,13
    condition = not (d==b or d==c or a==b or a==c or a==d or b==a+1 or c==a+1)
    result_assume = [(a, b, c, d, a), (a+1,b,c,d,a+1)]
    """
    config = {
        'setA': [7],
        'setB': list(range(1, 14)),
        'setC': list(range(1, 14)),
        'setD': [8],
        'which_mod': 1,
        'atypical_order': ['a', 'b', 'c', 'd', 'a'],
        'typical_order': ['a', 'b', 'c', 'd', 'a+1'],
        'result_assume': [
            ['a', 'b', 'c', 'd', 'a'],
            ['a+1', 'b', 'c', 'd', 'a+1']
        ],
        'condition': 'not (d==b or d==c or a==b or a==c or a==d or b==a+1 or c==a+1)'
    }
    
    def parse_value(value_str):
        """解析单个值，支持整数和有理数"""
        value_str = value_str.strip()
        
        # 如果是分数形式
        if '/' in value_str:
            try:
                numerator, denominator = map(int, value_str.split('/'))
                return numerator / denominator
            except:
                raise ValueError(f"无法解析分数: {value_str}")
        
        # 如果是小数形式
        elif '.' in value_str:
            try:
                return float(value_str)
            except:
                raise ValueError(f"无法解析小数: {value_str}")
        
        # 如果是整数
        else:
            try:
                return int(value_str)
            except:
                raise ValueError(f"无法解析整数: {value_str}")
    
    def parse_list(value_str):
        """解析列表，支持整数和有理数"""
        value_str = value_str.strip()
        
        # 移除可能的方括号
        if value_str.startswith('[') and value_str.endswith(']'):
            value_str = value_str[1:-1].strip()
        
        if not value_str:
            return []
        
        elements = []
        for item in value_str.split(','):
            item = item.strip()
            if item:
                elements.append(parse_value(item))
        
        return elements
    
    def parse_symbolic_list(value_str):
        """解析符号列表，如 [a, b, c, d, a]"""
        value_str = value_str.strip()
        
        if value_str.startswith('[') and value_str.endswith(']'):
            try:
                # 使用 ast.literal_eval 解析列表
                result = ast.literal_eval(value_str)
                if isinstance(result, list):
                    # 转换为字符串
                    return [str(item) for item in result]
            except:
                # 如果不能解析，尝试手动分割
                inner = value_str[1:-1].strip()
                if not inner:
                    return []
                return [x.strip() for x in inner.split(',') if x.strip()]
        
        # 逗号分隔格式
        return [x.strip() for x in value_str.split(',') if x.strip()]
    
    def parse_result_assume(value_str):
        """
        解析 result_assume 字段
        
        格式示例：
        result_assume = [(a, b, c, d, a), (a+1,b,c,d,a+1)]
        """
        value_str = value_str.strip()
        
        # 如果是列表格式
        if value_str.startswith('[') and value_str.endswith(']'):
            try:
                # 尝试使用 ast.literal_eval 解析
                # 注意：literal_eval 只能解析字面量，不能解析变量
                # 所以我们需要自定义解析
                
                # 使用正则表达式匹配元组
                # 移除外层方括号
                inner = value_str[1:-1].strip()
                
                result = []
                # 匹配格式如 (a, b, c, d, a) 或 (a+1,b,c,d,a+1)
                tuple_pattern = r'\(([^)]+)\)'
                tuples = re.findall(tuple_pattern, inner)
                
                for tup_str in tuples:
                    # 解析元组中的每个元素
                    elements = []
                    for elem in tup_str.split(','):
                        elem = elem.strip()
                        if elem:
                            elements.append(elem)
                    result.append(elements)
                
                return result
            except Exception as e:
                print(f"解析 result_assume 时出错: {e}")
                return config['result_assume']
        else:
            # 尝试直接解析
            try:
                # 简单的逗号分隔格式（不带元组）
                elements = []
                for item in value_str.split(','):
                    item = item.strip()
                    if item:
                        elements.append(item)
                return [elements]  # 包装成单个向量的列表
            except:
                return config['result_assume']
    
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                
                # 跳过注释和空行
                if not line or line.startswith('#'):
                    continue
                
                if '=' in line:
                    key, value = line.split('=', 1)
                    key = key.strip()
                    value = value.strip()
                    
                    try:
                        if key in ['setA', 'setB', 'setC', 'setD']:
                            config[key] = parse_list(value)
                            
                        elif key == 'which_mod':
                            config[key] = int(value)
                            
                        elif key in ['atypical_order', 'typical_order']:
                            # 解析符号列表
                            config[key] = parse_symbolic_list(value)
                            
                        elif key == 'result_assume':
                            config[key] = parse_result_assume(value)
                            
                        elif key == 'condition':
                            config[key] = value
                            
                        else:
                            print(f"第{line_num}行: 未知配置项 '{key}'，已忽略")
                            
                    except Exception as e:
                        print(f"第{line_num}行: 解析 '{key}' 时出错 - {e}")
                        print(f"  值: '{value}'")
        
        print(f"成功从文件 {filename} 读取配置")
        return config
        
    except FileNotFoundError:
        print(f"配置文件 {filename} 未找到，使用默认配置")
        input("程序暂停")
        return config
    except Exception as e:
        print(f"解析配置文件时出错: {e}")
        input("程序暂停")
        return config


def evaluate_condition(condition_str, variables):
    """
    安全地评估条件表达式
    variables: 包含变量值的字典，如 {'a': 1, 'b': 2, 'c': 3, 'd': 4}
    """
    try:
        # 使用ast进行基本的安全检查
        tree = ast.parse(condition_str, mode='eval')
        
        # 检查允许的变量名
        allowed_names = {'a', 'b', 'c', 'd', 'True', 'False', 'None'}
        
        for node in ast.walk(tree):
            if isinstance(node, ast.Name):
                if node.id not in allowed_names:
                    raise ValueError(f"不允许的变量名: {node.id}")
        
        # 安全地评估表达式
        result = eval(condition_str, {"__builtins__": {}}, variables)
        return bool(result)
        
    except Exception as e:
        print(f"评估条件时出错: {e}, 条件: {condition_str}")
        input("程序暂停,evaluate_condition")
        return False


def construct_vector(order, variables):
    """根据顺序构造向量"""
    values = []
    for param in order:
        if param == 'a':
            values.append(variables['a'])
        elif param == 'b':
            values.append(variables['b'])
        elif param == 'c':
            values.append(variables['c'])
        elif param == 'd':
            values.append(variables['d'])
        elif param == 'a+1':
            values.append(variables['a'] + 1)
        else:
            # 尝试解析表达式
            try:
                # 添加基本数学运算支持
                safe_locals = variables.copy()
                safe_locals.update({
                    'True': True,
                    'False': False,
                    'None': None
                })
                result = eval(param, {"__builtins__": {}}, safe_locals)
                values.append(result)
            except Exception as e:
                print(f"构造向量时出错: 参数 '{param}' - {e}")
                values.append(0)  # 默认值
                input("程序暂停,construct_vector")
    
    return values


def construct_assume_vectors(result_assume_spec, variables):
    """
    根据 result_assume 规范构造假设向量
    
    Args:
        result_assume_spec: result_assume 配置，如 [['a', 'b', 'c', 'd', 'a'], ['a+1', 'b', 'c', 'd', 'a+1']]
        variables: 变量字典，如 {'a': 1, 'b': 2, 'c': 3, 'd': 4}
    
    Returns:
        list of lists: 构造的向量列表
    """
    vectors = []
    
    for spec in result_assume_spec:
        values = []
        for expr in spec:
            try:
                # 安全地评估表达式
                safe_locals = variables.copy()
                safe_locals.update({
                    'True': True,
                    'False': False,
                    'None': None
                })
                
                # 评估表达式
                value = eval(expr, {"__builtins__": {}}, safe_locals)
                values.append(value)
            except Exception as e:
                print(f"构造假设向量时出错: 表达式 '{expr}' - {e}")
                values.append(0)  # 默认值
        
        vectors.append(values)
    
    return vectors


# 测试函数
def test_config_parser():
    """测试配置解析器"""
    
    # 创建测试配置文件
    config_content = """# 测试配置文件

# 集合定义
setA = [7]
setB = [1,2,3,4,5,6,7,8,9,10,11,12,13]
setC = [1,2,3,4,5,6,7,8,9,10,11,12,13]
setD = [8]

# 其他参数
which_mod = 1

# 向量顺序
atypical_order = [a, b, c, d, a]
typical_order = [a, b, c, d, a+1]

result_assume = [(a, b, c, d, a), (a+1,b,c,d,a+1)]

# 判断条件
condition = not (d==b or d==c or a==b or a==c or a==d or b==a+1 or c==a+1)
"""
    
    # 写入文件
    with open('test_config.txt', 'w', encoding='utf-8') as f:
        f.write(config_content)
    
    print("测试配置文件已创建: test_config.txt")
    print("=" * 60)
    
    # 解析配置
    config = parse_simple_config('test_config.txt')
    
    # 显示解析结果
    print("解析结果:")
    print(f"setA: {config['setA']}")
    print(f"setB: {config['setB'][:5]}...")
    print(f"setC: {config['setC'][:5]}...")
    print(f"setD: {config['setD']}")
    print(f"which_mod: {config['which_mod']}")
    print(f"atypical_order: {config['atypical_order']}")
    print(f"typical_order: {config['typical_order']}")
    print(f"result_assume: {config['result_assume']}")
    print(f"condition: {config['condition']}")
    
    print("\n" + "=" * 60)
    
    # 测试向量构造
    print("测试向量构造:")
    
    # 测试变量
    test_vars = {
        'a': 7,
        'b': 1,
        'c': 2,
        'd': 8
    }
    
    print(f"测试变量: {test_vars}")
    
    # 构造 atypical 向量
    atypical_vector = construct_vector(config['atypical_order'], test_vars)
    print(f"atypical_vector: {atypical_vector}")
    
    # 构造 typical 向量
    typical_vector = construct_vector(config['typical_order'], test_vars)
    print(f"typical_vector: {typical_vector}")
    
    # 构造 result_assume 向量
    assume_vectors = construct_assume_vectors(config['result_assume'], test_vars)
    print(f"\nresult_assume 向量:")
    for i, vec in enumerate(assume_vectors):
        print(f"  向量{i+1}: {vec}")
    
    print("\n" + "=" * 60)
    
    # 测试条件评估
    print("测试条件评估:")
    
    condition_result = evaluate_condition(config['condition'], test_vars)
    print(f"条件: {config['condition']}")
    print(f"评估结果: {condition_result}")
    
    # 测试几个不同的变量组合
    test_cases = [
        {'a': 7, 'b': 1, 'c': 2, 'd': 8, 'expected': True},
        {'a': 7, 'b': 7, 'c': 2, 'd': 8, 'expected': False},
        {'a': 6, 'b': 1, 'c': 7, 'd': 8, 'expected': False},
    ]
    
    print(f"\n条件测试:")
    for i, test in enumerate(test_cases, 1):
        result = evaluate_condition(config['condition'], test)
        status = "✓" if result == test['expected'] else "✗"
        print(f"{status} 测试 {i}: a={test['a']}, b={test['b']}, c={test['c']}, d={test['d']}")
        print(f"     结果: {result}, 期望: {test['expected']}")


if __name__ == "__main__":
    test_config_parser()
