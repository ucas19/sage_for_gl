from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.groups.abelian_gps.abelian_group import AbelianGroup
from sage.combinat.root_system.weyl_group import WeylGroup
from sage.matrix.constructor import diagonal_matrix, identity_matrix
from sage.modules.free_module_element import vector
from SemidirectProductGroup import SemidirectProductGroup
from sage.combinat.root_system.cartan_type import CartanType
#from sage.groups.artin import CoxeterGroup
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.all import CoxeterGroup, AA
from collections import Counter
from sage.all import zero_vector

def check_arrays(A, B, n, m):
    def remove_matches(front, back):
        # front, back 是列表
        front_ctr = Counter(front)
        back_ctr = Counter(back)
        removed_count = 0
        
        # 遍历前段的元素（按值，但考虑重复）
        # 为了避免修改迭代中的字典，我们复制键的列表
        for x in list(front_ctr.keys()):
            abs_x = abs(x)
            # 在后段找绝对值相同的元素
            # 需要找后段中绝对值等于abs_x的任意一个元素（比如y，|y|==abs_x）
            # 但可能有多个不同符号的，我们只关心绝对值
            # 所以需要遍历后段键找匹配绝对值的
            match_found = None
            for y in list(back_ctr.keys()):
                if abs(y) == abs_x:
                    match_found = y
                    break
            if match_found is not None:
                # 能匹配，删掉一对
                min_pair = min(front_ctr[x], back_ctr[match_found])
                front_ctr[x] -= min_pair
                back_ctr[match_found] -= min_pair
                removed_count += 2 * min_pair
                if front_ctr[x] == 0:
                    del front_ctr[x]
                if back_ctr[match_found] == 0:
                    del back_ctr[match_found]
        return front_ctr, back_ctr, removed_count

    # 对A处理
    front_A = A[:n]
    back_A = A[n:]
    front_A_ctr_rem, back_A_ctr_rem, removed_A = remove_matches(front_A, back_A)
    
    # 对B处理
    front_B = B[:n]
    back_B = B[n:]
    front_B_ctr_rem, back_B_ctr_rem, removed_B = remove_matches(front_B, back_B)
    
    if removed_A != removed_B:
        return False
    
    # 比较剩余部分（绝对值多重集合）
    def abs_multiset(ctr):
        # 将计数器中的键取绝对值，并展开成多重集合列表（其实用Counter合并绝对值相同的键）
        abs_ctr = Counter()
        for k, v in ctr.items():
            abs_ctr[abs(k)] += v
        return abs_ctr
    
    abs_front_A = abs_multiset(front_A_ctr_rem)
    abs_back_A = abs_multiset(back_A_ctr_rem)
    abs_front_B = abs_multiset(front_B_ctr_rem)
    abs_back_B = abs_multiset(back_B_ctr_rem)
    
    if abs_front_A == abs_front_B and abs_back_A == abs_back_B:
        return True
    return False

def sym_3( n, m):

    V_0 = []
    V_1 = []
    result = []
    v = zero_vector(QQ, n+m)
    V_0.append(v[:])
    for i in range(m):
        v = zero_vector(QQ, n+m)
        v[n+i] = 1
        V_0.append(v[:])
        v[n+i] = -1
        V_0.append(v[:])

    for i in range(n):
        v = zero_vector(QQ, n+m)
        v[i] = 1
        V_1.append(v[:])
        v[i] = -1
        V_1.append(v[:])
    
    for i in range(len(V_1)):
        for j in range(i+1,len(V_1)):
            for k in range(j+1,len(V_1)):
                tem = V_1[i]+V_1[j]+V_1[k]
                result.append(tem[:])

    for i in range(len(V_1)):
        for j in range(i+1,len(V_1)):
            for k in range(len(V_0)):
                tem = V_1[i]+V_1[j]+V_0[k]
                result.append(tem[:])

    for i in range(len(V_0)):
        for j in range(i,len(V_0)):
            for k in range(len(V_1)):
                tem = V_0[i]+V_0[j]+V_1[k]
                result.append(tem[:])

    for i in range(len(V_0)):
        for j in range(i,len(V_0)):
            for k in range(j,len(V_0)):
                tem = V_0[i]+V_0[j]+V_0[k]
                result.append(tem[:])
    return result

def wedge_3( n, m):

    V_0 = []
    V_1 = []
    result = []
    v = zero_vector(QQ, n+m)
    V_0.append(v[:])
    for i in range(m):
        v = zero_vector(QQ, n+m)
        v[n+i] = 1
        V_0.append(v[:])
        v[n+i] = -1
        V_0.append(v[:])

    for i in range(n):
        v = zero_vector(QQ, n+m)
        v[i] = 1
        V_1.append(v[:])
        v[i] = -1
        V_1.append(v[:])
    
    for i in range(len(V_0)):
        for j in range(i+1,len(V_0)):
            for k in range(j+1,len(V_0)):
                tem = V_0[i]+V_0[j]+V_0[k]
                result.append(tem[:])

    for i in range(len(V_0)):
        for j in range(i+1,len(V_0)):
            for k in range(len(V_1)):
                tem = V_0[i]+V_0[j]+V_1[k]
                result.append(tem[:])


    for i in range(len(V_1)):
        for j in range(i,len(V_1)):
            for k in range(len(V_0)):
                tem = V_1[i]+V_1[j]+V_0[k]
                result.append(tem[:])

    for i in range(len(V_1)):
        for j in range(i,len(V_1)):
            for k in range(j,len(V_1)):
                tem = V_1[i]+V_1[j]+V_1[k]
                result.append(tem[:])
    return result





class Lowest_Module:
    
    def __init__(self, n, m):
        self.V = []
        self.S2V = []
        self.g = []
        self.S3V = []
        self.W3V = wedge_3(n,m)
        self.S3VV = sym_3(n,m)
        self.basis_plus = []
        
        for i in range(n+m-1):
            v = zero_vector(QQ, n+m)
            v[i] = 1
            v[i+1] = -1
            self.basis_plus.append(v[:])
        v = zero_vector(QQ, n+m)
        v[n+m-1] = 1
        self.basis_plus.append(v[:])



        v = zero_vector(QQ, n+m)
        self.V.append(v[:])

        for i in range(n):
            v = zero_vector(QQ, n+m)
            v[i] = 1
            self.V.append(v[:])

            v[i] = -1
            self.V.append(v[:])

        for i in range(m):
            v = zero_vector(QQ, n+m)
            v[n+i] = 1
            self.V.append(v[:])
            v[n+i] = -1
            self.V.append(v[:])

        v = zero_vector(QQ, n+m)
        self.S2V.append(v[:])
        for i in range(n+m):
            self.g.append(v[:])
            self.S2V.append(v[:])
        for i in range(n):
            v = zero_vector(QQ, n+m)
            v[i] = 2
            self.g.append(v[:])
            v[i] = -2
            self.g.append(v[:])
            v[i] = 1
            self.g.append(v[:])
            self.S2V.append(v[:])
            v[i] = -1
            self.g.append(v[:])
            self.S2V.append(v[:])
        for i in range(n):
            for j in range(i+1,n):
                v = zero_vector(QQ, n+m) 
                v[i] = 1
                v[j] = -1
                self.g.append(v[:])
                self.S2V.append(v[:])
                v[i] = -1
                v[j] = 1
                self.g.append(v[:])
                self.S2V.append(v[:])
                v[i] = 1
                v[j] = 1
                self.g.append(v[:])
                self.S2V.append(v[:])
                v[i] = -1
                v[j] = -1
                self.g.append(v[:])
                self.S2V.append(v[:])



        for i in range(m):
            v = zero_vector(QQ, n+m)
            v[n+i] = 1
            self.g.append(v[:])
            self.S2V.append(v[:])
            v[n+i] = -1
            self.g.append(v[:])
            self.S2V.append(v[:])
            v[n+i] = 2
            self.S2V.append(v[:])
            v[n+i] = -2
            self.S2V.append(v[:])
        for i in range(m):
            for j in range(i+1,m):
                v = zero_vector(QQ,n+m)
                v[n+i] = 1
                v[n+j] = -1
                self.g.append(v[:])
                self.S2V.append(v[:])
                v[n+i] = -1
                v[n+j] = 1
                self.g.append(v[:])
                self.S2V.append(v[:])
                v[n+i] = 1
                v[n+j] = 1
                self.g.append(v[:])
                self.S2V.append(v[:])
                v[n+i] = -1
                v[n+j] = -1
                self.g.append(v[:])
                self.S2V.append(v[:])
        for i in range(n):
            for j in range(m):
                v = zero_vector(QQ,n+m)
                v[i] = 1
                v[n+j] = -1
                self.g.append(v[:])
                self.S2V.append(v[:])
                v[i] = -1
                v[n+j] = 1
                self.g.append(v[:])
                self.S2V.append(v[:])
                v[i] = 1
                v[n+j] = 1
                self.g.append(v[:])
                self.S2V.append(v[:])
                v[i] = -1
                v[n+j] = -1
                self.g.append(v[:])
                self.S2V.append(v[:])
        



        for i in range(m+1+n):
            v = zero_vector(QQ,n+m)
            self.S3V.append(v[:])

        for i in range(m):
            for j in range(m+n+1):
                v = zero_vector(QQ,n+m)
                v[n+i] = 1
                self.S3V.append(v[:])
                v[n+i] = -1
                self.S3V.append(v[:])

        for i in range(m):
            v = zero_vector(QQ,n+m)
            v[n+i] = 2
            self.S3V.append(v[:])
            v[n+i] = -2
            self.S3V.append(v[:])
            v[n+i] = 3
            self.S3V.append(v[:])
            v[n+i] = -3
            self.S3V.append(v[:])

        for i in range(m):
            for j in range(i+1,m):
                v = zero_vector(QQ,n+m)
                v[n+i] = 1
                v[n+j] = -1
                self.S3V.append(v[:])
                v[n+i] = -1
                v[n+j] = 1
                self.S3V.append(v[:])
                v[n+i] = 1
                v[n+j] = 1
                self.S3V.append(v[:])
                v[n+i] = -1
                v[n+j] = -1
                self.S3V.append(v[:])
                
        for i in range(m):
            for j in range(i+1,m):
                for k in range(j+1,m):
                    syms = [1,-1]
                    for sym1 in syms:
                        for sym2 in syms:
                            for sym3 in syms:
                                v = zero_vector(QQ,n+m)
                                v[n+i] = sym1
                                v[n+j] = sym2
                                v[n+k] = sym3
                                self.S3V.append(v[:])
        for i in range(m):
            for j in range(m):
                if i==j:
                    continue
                syms = [1,-1]
                for sym1 in syms:
                    for sym2 in syms:
                        v = zero_vector(QQ,n+m)
                        v[n+i] = 2 *sym1
                        v[n+j] = sym2
                        self.S3V.append(v[:])

        for i in range(n):
            for j in range(m+n):
                v = zero_vector(QQ,n+m)
                v[i] = 1
                self.S3V.append(v[:])
                v[i] = -1
                self.S3V.append(v[:])

        for i in range(n):
            for j in range(i+1,n):
                syms = [1,-1]
                for sym1 in syms:
                    for sym2 in syms:
                        v = zero_vector(QQ,n+m)
                        v[i] = sym1
                        v[j] = sym2
                        self.S3V.append(v[:])

        for i in range(n):
            for j in range(i+1,n):
                for k in range(j+1,n):
                    syms = [1,-1]
                    for sym1 in syms:
                        for sym2 in syms:
                            for sym3 in syms:
                                v = zero_vector(QQ,n+m)
                                v[i] = sym1
                                v[j] = sym2
                                v[k] = sym3
                                self.S3V.append(v[:])

        for i in range(n):
            for j in range(m):
                syms = [1,-1]
                for sym1 in syms:
                    for sym2 in syms:
                        v = zero_vector(QQ,n+m)
                        v[i] = sym1
                        v[n+j] = sym2
                        self.S3V.append(v[:])
                        v[i] = sym1
                        v[n+j] = 2*sym2
                        self.S3V.append(v[:])

        for i in range(n):
            for j in range(i+1,n):
                for k in range(m):
                    syms = [1,-1]
                    for sym1 in syms:
                        for sym2 in syms:
                            for sym3 in syms:
                                v = zero_vector(QQ,n+m)
                                v[i] = sym1
                                v[j] = sym2
                                v[n+k] = sym3
                                self.S3V.append(v[:])
        for i in range(m):
            for j in range(i+1,m):
                for k in range(n):
                    syms = [1,-1]
                    for sym1 in syms:
                        for sym2 in syms:
                            for sym3 in syms:
                                v = zero_vector(QQ,n+m)
                                v[n+i] = sym1
                                v[n+j] = sym2
                                v[k] = sym3
                                self.S3V.append(v[:])




        
