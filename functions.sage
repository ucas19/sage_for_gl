from typing import Counter
from Mu_with_message import Mu_With_Message
from weyl_group_Bn import Weyl_Group_Bn
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.rings.rational_field import QQ
from lowest_module import Lowest_Module
from integer_com import is_nonnegative_integer_combination_simple
from sage_integer_com import is_nonnegative_integer_combination_sage
from Mu_with_message import *
from collections import defaultdict
from sage.all import *
    

def get_coset_representatives(W_Group, J):
    simple_refls = W_Group.W.simple_reflections()
    gens = []
    gens = [simple_refls[i] for i in J]
    parabolic_subgroup = W_Group.W.subgroup(gens)
    #    print(parabolic_subgroup)
    #print("Order of parabolic subgroup:", parabolic_subgroup.order())
    min_left_coset_reps = [w for w in W_Group.W if all(not w.has_right_descent(i) for i in J)]
    #print("Minimal left coset representatives:")
    #for rep in min_left_coset_reps:
    #    print(rep.reduced_word(), "length:", rep.length())
    #print("Number of minimal left coset representatives:", len(min_left_coset_reps))
    return min_left_coset_reps

def get_S(Weyl_W, lambda_sp_next):
    #找到抛物子群代表元
    print(f"get_S({lambda_sp_next})")

    S_sp = [] # 这是W_sp的子抛物子群生成元集合（序号）
#    for w in W_sp.W.gens():
#        print(f"{w} ---- {w.to_matrix()}")
#        if w.to_matrix() * lambda_sp_next == lambda_sp_next:
    #            print(f"{w} -- {w.to_matrix()}")

    #这里是weyl的所有生成元集
    S_of_W_sp = Weyl_W.simple_reflections()
    #    print(S_of_W_sp[i].coset_representative([2],side = "left"))
    for i in S_of_W_sp.keys():
        #        print(f"{i}-- --{S_of_W_sp[i].to_matrix()}")
        if S_of_W_sp[i].to_matrix() * lambda_sp_next == lambda_sp_next :
            S_sp.append(i)
    return S_sp



def calc_w_mu( W_Group , lambda_before ):
    n = W_Group.n
    w_result = []
    lambda_result = lambda_before
    calc_sum = 0
    for w in W_Group.W:
        flag = 0
        mu = w.to_matrix() * lambda_before
 #       print(f"{mu} = {lambda_before}*\n {w.to_matrix()}")
        for i in range(n):
            for j in range(i+1,n):
                delta_add_delta = [0 for _ in range(n)]
                delta_add_delta[i] = 1
                delta_add_delta[j] = -1
                sum = 0
                for k in range(n):
                    sum = sum + delta_add_delta[k] * mu[k]
 #                   print(f"e{delta_add_delta} * {mu} = {sum}")
                if sum > 0:
 #                   print(f"b{delta_add_delta} * {mu} = {sum}")
                    flag = 1
                    break

        if flag == 0:
            lambda_result = w.to_matrix() * lambda_before
            w_result.append(w)
  #          print(f"{w} and {lambda_result}")
            calc_sum = calc_sum+1
    #    print(calc_sum)
    if calc_sum == 1:
        #说明lambda_before正则
        return w_result[0].inverse(), lambda_result
    elif calc_sum > 1 :
        #这种情况下
        #说明lambda_before正则
        S_set = get_S(W_Group.W, lambda_result) 
        for w in w_result:
            w_inverse = w.inverse()
            if all(not w_inverse.has_right_descent(i) for i in S_set):
                 w_result_only = w_inverse
        return w_result_only ,lambda_result

    else:
        print(f"{calc_sum},{w_result}")
        print("error计算错误")
        return None
def minimal_coset_representatives_manual(W, S):
    """
    手动计算最小陪集代表元
    """
    # 创建抛物子群
    W_S = W.parabolic_subgroup(S)
    # 获取所有左陪集
    cosets = W.cosets(W_S, side='left')
    # 对每个陪集找到长度最小的元素
    min_reps = []
    for coset in cosets:
        # 按长度排序并取第一个（最小长度）
        min_rep = min(coset, key=lambda w: w.length())
        min_reps.append(min_rep)
    
    return min_reps



def remove_matches(front, back):
    # front, back 是列表
    removed_count = 0
    front_tem = front[:]
    for x in front_tem:
        if x in back:
            front.remove(x)
            back.remove(x)
            removed_count += 1
    front_ctr = Counter(front) 
    back_ctr = Counter(back)
    return front[:],back[:], front_ctr, back_ctr, removed_count
    """
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
    """
def abs_multiset(ctr):
    # 将计数器中的键取绝对值，并展开成多重集合列表（其实用Counter合并绝对值相同的键）
    abs_ctr = Counter()
    for k, v in ctr.items():
        abs_ctr[abs(k)] += v
    return abs_ctr

def check_arrays(A, B, n, m):

    # 对A处理
    front_A = A[:n]
    back_A = A[n:]
    f_A,b_A, front_A_ctr_rem, back_A_ctr_rem, removed_A = remove_matches(front_A, back_A)
    
    # 对B处理
    front_B = B[:n]
    back_B = B[n:]
    f_B,b_B,front_B_ctr_rem, back_B_ctr_rem, removed_B = remove_matches(front_B, back_B)
    
    if removed_A != removed_B:
        return False
    
    # 比较剩余部分（绝对值多重集合）
    
    abs_front_A = front_A_ctr_rem
    abs_back_A = back_A_ctr_rem
    abs_front_B = front_B_ctr_rem
    abs_back_B = back_B_ctr_rem
    
    if abs_front_A == abs_front_B and abs_back_A == abs_back_B:
        return True
    return False
    

def P_tensor_V(lambda_sp_plus_so,sum_sp_plus_so,V,n,m):

    print("P_mu_tensor_V如下:")
    P_mu_tensor_V_befor_Pr = []
    P_mu_tensor_V_after_Pr = []
    calc_sum = 0
    for v in sum_sp_plus_so:
        for w in V:
            P_mu_tensor_V_befor_Pr.append(v+w)
            #print(f"{calc_sum}:  {v+w}")
            calc_sum +=1
            if(check_arrays((v+w).list(), lambda_sp_plus_so.list(), n, m)):
                P_mu_tensor_V_after_Pr.append(v+w)
    print(f"不再展示P_mu_tensor_V了,个数是:{calc_sum}")
    return P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr



def P_tensor_V_not_show(lambda_sp_plus_so,sum_sp_plus_so,V,n,m):

#    print("P_mu_tensor_V如下:")
    P_mu_tensor_V_befor_Pr = []
    P_mu_tensor_V_after_Pr = []
    calc_sum = 1
    for v in sum_sp_plus_so:
        for w in V:
            P_mu_tensor_V_befor_Pr.append(v+w)
            #print(f"{calc_sum}:  {v+w}")
            calc_sum +=1
            if(check_arrays((v+w).list(), lambda_sp_plus_so.list(), n, m)):
                P_mu_tensor_V_after_Pr.append(v+w)
#    print(f"不再展示P_mu_tensor_V了,个数是:{calc_sum}")
    return P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr




def which_one_lowest(P_mu_tensor_V_after_Pr, basis_plus):

    
    vector_sets = P_mu_tensor_V_after_Pr[:]
    show_P_mu = []
    not_result = []

    for lambda_sp_plus_so in vector_sets:
        if lambda_sp_plus_so in not_result:
            continue
        flag = True
        for v in vector_sets:   
            result, coefficients = is_nonnegative_integer_combination_sage(lambda_sp_plus_so-v, basis_plus)
            if result:
                flag = False
#                print(f"测试1-{lambda_sp_plus_so}")
#                print(f"测试2-{v}")
#                print(f"测试3-{coefficients}")
                break
            
            result, coefficients = is_nonnegative_integer_combination_sage(v-lambda_sp_plus_so, basis_plus)
            if result:
                if v not in not_result:
                    not_result.append(v)

        if flag:   
            show_P_mu.append(lambda_sp_plus_so)
    return show_P_mu

def equiva_vec(v, n, m):
    A = v.list()
    front_A = A[:n]
    back_A = A[n:]
    front,back,front_A_ctr_rem, back_A_ctr_rem, removed_A = remove_matches(front_A, back_A)
    part1=tuple(sorted(front))
    part2=tuple(sorted(back))
    return (part1,part2)

def selete_block(P_mu_tensor_V_befor_Pr,n,m):

    groups = defaultdict(list)
    for v in P_mu_tensor_V_befor_Pr:
        key = equiva_vec(v,n,m)
        #print(f"{calc_sum}:  {v+w}")
        groups[key].append(v)
    return groups            


def P_tensor_V_show(sum_sp_plus_so,which_mod,n,m):
    
    lowest_module = Lowest_Module(n,m)
    V = lowest_module.get_module(which_mod)

    print("P_mu_tensor_V如下:")
    P_mu_tensor_V_befor_Pr = []
    calc_sum = 0
    for v in sum_sp_plus_so:
        for w in V:
            P_mu_tensor_V_befor_Pr.append(v+w)
            #print(f"{calc_sum}:  {v+w}")
            calc_sum +=1
    
    print(f"不再展示P_mu_tensor_V了,个数是:{calc_sum}")

    select = selete_block(P_mu_tensor_V_befor_Pr,n,m)
    groups_num = 1
    for k, group in select.items():
        print(f"--------第{groups_num}组---------")
        groups_num += 1
        lowest_weights = which_one_lowest(group,lowest_module.basis_plus)
        for i in range(len(lowest_weights)):
            print(f"{i+1}: {lowest_weights[i]}")
    return select


def is_tensor_V_true(lambda_sp_plus_so, P_mu_tensor_V_after_Pr, basis_plus):
    print("Pr_P_mu_tensor_V如下:")
    calc_sum = 1
    flag = True
    is_in = False
    immutable_vecs = [vector(v, immutable=True) for v in P_mu_tensor_V_after_Pr]
    P_mu_tensor_V_after_Pr_dic = Counter(immutable_vecs)
    show_P_mu = []
    for v, k in P_mu_tensor_V_after_Pr_dic.items():
        if v == lambda_sp_plus_so:
            is_in = True
        print(f"{calc_sum}:  {v} 数量:{k}")
        calc_sum +=1
        result, coefficients = is_nonnegative_integer_combination_sage(lambda_sp_plus_so-v, basis_plus)
        if result:
            print(f"第{calc_sum-1}个向量{v}比{lambda_sp_plus_so}小, 因此空间不成立系数是{coefficients}")
            flag = False
        else:
            show_P_mu.append(v)
    if is_in== False:
        print(f"数据有问题，{lambda_sp_plus_so}不在里面。")
        return False
    if flag:
#        print(f"整理如下:")
#        print(show_P_mu)
        print(f"展开整理如下:")
        print(P_mu_tensor_V_after_Pr)
        print(f"数量: {len(P_mu_tensor_V_after_Pr)}")

    return flag


def is_tensor_V_true_not_show(lambda_sp_plus_so, P_mu_tensor_V_after_Pr, basis_plus):
#    print("Pr_P_mu_tensor_V如下:")
    calc_sum = 1
    flag = True
    is_in = False
    immutable_vecs = [vector(v, immutable=True) for v in P_mu_tensor_V_after_Pr]
    P_mu_tensor_V_after_Pr_dic = Counter(immutable_vecs)
    show_P_mu = []
    for v, k in P_mu_tensor_V_after_Pr_dic.items():
        if v == lambda_sp_plus_so:
            is_in = True
#        print(f"{calc_sum}:  {v} 数量:{k}")
        calc_sum +=1
        
        result, coefficients = is_nonnegative_integer_combination_sage(lambda_sp_plus_so-v, basis_plus)

#        test = vector(QQ,[-1/2,-3/2,3/2,3/2])
#        if test == v:
#            print()

        if result:
 #           print(f"第{calc_sum-1}个向量{v}比{lambda_sp_plus_so}小, 因此空间不成立系数是{coefficients}")
            flag = False
        else:
            show_P_mu.append(v)
    if is_in== False:
 #       print(f"数据有问题，{lambda_sp_plus_so}不在里面。")
        return False
#    if flag:
 #       print(f"整理如下:")
  #      print(show_P_mu)
#        print(f"展开整理如下:")
#        print(P_mu_tensor_V_after_Pr)
#        print(f"数量: {len(P_mu_tensor_V_after_Pr)}")

    return flag



def test_kl(n):


    W_sp = Weyl_Group_Bn(n)

    R.<q> = LaurentPolynomialRing(QQ)
    KL_sp = KazhdanLusztigPolynomial(W_sp.W, q)
    
    print("111111")
    print("验证是否都是1")
    for rep1 in W_sp.W:
        for rep2 in W_sp.W:
        #        print(f"这里是rep--{rep}")
            if rep1.bruhat_le(rep2):

                k_l_poly = KL_sp.P(rep1,rep2)
                k_l_on_one=k_l_poly(1)
                if not k_l_on_one == 1:
                    print(f"不是1:{rep1},{rep2}")
    print("222222")
    print("验证是否都是0")
    for rep1 in W_sp.W:
        for rep2 in W_sp.W:
            if rep2.bruhat_le(rep1) and not rep1 == rep2:

                k_l_poly = KL_sp.P(rep1,rep2)
                k_l_on_one=k_l_poly(1)
                if not k_l_on_one == 0:
                    print(f"不是2:{rep1},{rep2}={k_l_on_one}")


def K_L_decompose_no_kl(W_sp,w_sp, lambda_sp_next, W_so, w_so, lambda_so_next):

    S_sp = [] # 这是W_sp的子抛物子群生成元集合（序号）
    S_so = [] # 这是W_so的子抛物子群生成元集合（序号）
#    for w in W_sp.W.gens():
#        print(f"{w} ---- {w.to_matrix()}")
#        if w.to_matrix() * lambda_sp_next == lambda_sp_next:
    #            print(f"{w} -- {w.to_matrix()}")

    #这里是weyl的所有生成元集
    S_of_W_sp = W_sp.W.simple_reflections()
    S_of_W_so = W_so.W.simple_reflections()
    #    print(S_of_W_sp[i].coset_representative([2],side = "left"))
    for i in S_of_W_sp.keys():
        #        print(f"{i}-- --{S_of_W_sp[i].to_matrix()}")
        if S_of_W_sp[i].to_matrix() * lambda_sp_next == lambda_sp_next :
            S_sp.append(i)

    for i in S_of_W_so.keys():
        if S_of_W_so[i].to_matrix() * lambda_so_next == lambda_so_next :
            S_so.append(i)

    #    for i in range(W_sp.W.rank()):
    #        print(f"{W_sp.W[i].to_matrix()} * {lambda_sp_next}")
    #        if W_sp.W[i].to_matrix() * lambda_sp_next == lambda_sp_next :
    #            S_sp.append(i+1)
    # print(f"{W_sp.W[i].to_matrix()} * {lambda_sp_next}")
    #    print(S_sp)
    #    min_reps_sp = W_sp.W_cox.coset_representative(S_sp,side = "left" ) 
    #    W_S = W_sp.W.parabolic_subgroup(S_of_W_sp)
    #这里得到S_sp（S_so）为生成元的子抛物子群的最小左陪集代表元集合
    min_left_coset_reps_sp = get_coset_representatives(W_sp, S_sp)
    min_left_coset_reps_so = get_coset_representatives(W_so, S_so)

    #这里是在最小左陪集代表元集合中，找到符合条件的，也就是比w_sp大的元素
    sum_sp_weyl = []
    sum_so_weyl = []
    

    #W_sp_coxeter = W_sp.W_cox
    #W_so_coxeter = W_so.W_cox


    for rep in min_left_coset_reps_sp:
        #        print(f"这里是rep--{rep}")
        is_bruhat_le_reverse = w_sp.bruhat_le(rep)
        if w_sp.bruhat_le(rep) or w_sp == rep:
            #            print(f"sp的kl多项式的sum中weyl元={rep}")
            #print(f"sp的kl多项式的sum中后M_?集合={rep.to_matrix() * lambda_sp_next }")
            #           sum_sp.append(rep.to_matrix() * lambda_sp_next)

            #w_sp_cox = W_sp_coxeter.from_reduced_word(w_sp.reduced_word())
            #rep_cox = W_sp_coxeter.from_reduced_word(rep.reduced_word())
            #k_l_on_one = W_sp_coxeter.kazhdan_lusztig_polynomial(rep_cox*W_sp_coxeter.long_element(),w_sp_cox*W_sp_coxeter.long_element())
            sum_sp_weyl.append(rep)

    for rep in min_left_coset_reps_so:
        #        print(f"这里是rep--{rep}")
        is_bruhat_le_reverse = w_so.bruhat_le(rep)
        if w_so.bruhat_le(rep):
            #print(f"so的kl多项式的sum中weyl元={rep}")
            #print(f"so的kl多项式的sum中后M_?集合={rep.to_matrix() * lambda_so_next }")
            #           sum_so.append(rep.to_matrix() * lambda_so_next)
            #w_so_cox = W_so_coxeter.from_reduced_word(w_so.reduced_word())
            #rep_cox = W_so_coxeter.from_reduced_word(rep.reduced_word())
            #k_l_on_one = W_so_coxeter.kazhdan_lusztig_polynomial(rep_cox*W_so_coxeter.long_element(),w_so_cox*W_so_coxeter.long_element())
            sum_so_weyl.append(rep)

    #将sp和so连起来
    sum_sp_plus_so = []
    #    print("P_mu的分解如下")
    calc_sum = 1
    for v in sum_sp_weyl:
        for w in sum_so_weyl:
            v_plus_w = vector(QQ, list(v.to_matrix()*lambda_sp_next)+list(w.to_matrix()*(-lambda_so_next)))
            sum_sp_plus_so.append(v_plus_w)
            #            print(f"{calc_sum}: {v_plus_w}")
            calc_sum += 1
    return sum_sp_weyl, sum_so_weyl, sum_sp_plus_so


def K_L_decompose(W_sp,w_sp, lambda_sp_next, W_so, w_so, lambda_so_next):

    S_sp = [] # 这是W_sp的子抛物子群生成元集合（序号）
    S_so = [] # 这是W_so的子抛物子群生成元集合（序号）
#    for w in W_sp.W.gens():
#        print(f"{w} ---- {w.to_matrix()}")
#        if w.to_matrix() * lambda_sp_next == lambda_sp_next:
    #            print(f"{w} -- {w.to_matrix()}")

    #这里是weyl的所有生成元集
    S_of_W_sp = W_sp.W.simple_reflections()
    S_of_W_so = W_so.W.simple_reflections()
    #    print(S_of_W_sp[i].coset_representative([2],side = "left"))
    for i in S_of_W_sp.keys():
        #        print(f"{i}-- --{S_of_W_sp[i].to_matrix()}")
        if S_of_W_sp[i].to_matrix() * lambda_sp_next == lambda_sp_next :
            S_sp.append(i)

    for i in S_of_W_so.keys():
        if S_of_W_so[i].to_matrix() * lambda_so_next == lambda_so_next :
            S_so.append(i)
    #    for i in range(W_sp.W.rank()):
    #        print(f"{W_sp.W[i].to_matrix()} * {lambda_sp_next}")
    #        if W_sp.W[i].to_matrix() * lambda_sp_next == lambda_sp_next :
    #            S_sp.append(i+1)
    # print(f"{W_sp.W[i].to_matrix()} * {lambda_sp_next}")
    #    print(S_sp)
    #    min_reps_sp = W_sp.W_cox.coset_representative(S_sp,side = "left" ) 
    #    W_S = W_sp.W.parabolic_subgroup(S_of_W_sp)
    #这里得到S_sp（S_so）为生成元的子抛物子群的最小左陪集代表元集合
    min_left_coset_reps_sp = get_coset_representatives(W_sp, S_sp)
    min_left_coset_reps_so = get_coset_representatives(W_so, S_so)
    print("S_sp")
    print("min_left_coset_reps_sp")
    print("S_so")
    print("min_left_coset_reps_so")

    #这里是在最小左陪集代表元集合中，找到符合条件的，也就是比w_sp大的元素
    sum_sp_weyl = []
    sum_so_weyl = []
    

    #W_sp_coxeter = W_sp.W_cox
    #W_so_coxeter = W_so.W_cox

    R.<q> = LaurentPolynomialRing(QQ)
    KL_sp = KazhdanLusztigPolynomial(W_sp.W, q)
 #   KL.P(s2,s3*s2*s3*s1*s2) 
    
    KL_so = KazhdanLusztigPolynomial(W_so.W, q)

 
   #这里是是测试
    #    elements = list(W_sp.W)
    #    n = len(elements)
    #    print("测试w")
    #    print(elements)
    #    kl_matrix = matrix(ZZ, n, n)  # KL多项式在1处的值通常是整
    #
    #    for i, v in enumerate(elements):
    #        for j, w in enumerate(elements):
    #            if v.bruhat_le(w):
    #                result = KL_sp.P(v, w)
    #                kl_matrix[i, j] = result(1)
    #    print("KL多项式在 q=1 处的取值矩阵:")
    #    print(kl_matrix)

 

    for rep in min_left_coset_reps_sp:
        #        print(f"这里是rep--{rep}")
        is_bruhat_le_reverse = w_sp.bruhat_le(rep)
        if w_sp.bruhat_le(rep) or w_sp == rep:
            #sum_sp.append(rep.to_matrix() * lambda_sp_next)

            #w_sp_cox = W_sp_coxeter.from_reduced_word(w_sp.reduced_word())
            #rep_cox = W_sp_coxeter.from_reduced_word(rep.reduced_word())
            #k_l_on_one = W_sp_coxeter.kazhdan_lusztig_polynomial(rep_cox*W_sp_coxeter.long_element(),w_sp_cox*W_sp_coxeter.long_element())
            k_l_poly = KL_sp.P(rep*W_sp.W.long_element(),w_sp*W_sp.W.long_element())
            k_l_on_one=k_l_poly(1)
#            if k_l_on_one > 1:

            print(f"sp: {w_sp}")
            print(f"sp_kl: {k_l_poly}")
            print(f"sp的kl多项式的sum中weyl元={rep}")
            print(f"sp的kl多项式的sum中后M_?集合={rep.to_matrix() * lambda_sp_next }")

            for i in range(k_l_on_one):
                sum_sp_weyl.append(rep)

    for rep in min_left_coset_reps_so:
        #        print(f"这里是rep--{rep}")
        is_bruhat_le_reverse = w_so.bruhat_le(rep)
        if w_so.bruhat_le(rep):
            #print(f"so的kl多项式的sum中weyl元={rep}")
            #print(f"so的kl多项式的sum中后M_?集合={rep.to_matrix() * lambda_so_next }")
            #           sum_so.append(rep.to_matrix() * lambda_so_next)
            #w_so_cox = W_so_coxeter.from_reduced_word(w_so.reduced_word())
            #rep_cox = W_so_coxeter.from_reduced_word(rep.reduced_word())
            #k_l_on_one = W_so_coxeter.kazhdan_lusztig_polynomial(rep_cox*W_so_coxeter.long_element(),w_so_cox*W_so_coxeter.long_element())
            k_l_poly = KL_so.P(rep*W_so.W.long_element(),w_so*W_so.W.long_element())
            k_l_on_one= k_l_poly(1)

            print(f"so: {w_so}")
            print(f"so_kl: {k_l_poly}")
            print(f"so的kl多项式的sum中weyl元={rep}")
            print(f"so的kl多项式的sum中后M_?集合={rep.to_matrix() * lambda_so_next }")
            for i in range(k_l_on_one):
                sum_so_weyl.append(rep)

    #将sp和so连起来
    sum_sp_plus_so = []
    #    print("P_mu的分解如下")
    calc_sum = 1
    for v in sum_sp_weyl:
        for w in sum_so_weyl:
            v_plus_w = vector(QQ, list(v.to_matrix()*lambda_sp_next)+list(w.to_matrix()*(-lambda_so_next)))
            sum_sp_plus_so.append(v_plus_w)
            #        print(f"{calc_sum}: {v_plus_w}")
            calc_sum += 1
    return sum_sp_weyl, sum_so_weyl, sum_sp_plus_so




def judge_mu_in_P(lambda_init,n,m,depth):
#由于算力问题，只能重复五次，应该是够用的

    #这里是所有偶正根
    Phi_0_plus = []
    for i in range(n):
        for j in range(i+1,n):
            v = zero_vector(QQ,n+m)
            v[i] = 1
            v[j] = -1
            Phi_0_plus.append(v[:])

    for i in range(m):
        for j in range(i+1,m):
            v = zero_vector(QQ,n+m)
            v[n+i] = -1
            v[n+j] = 1
            Phi_0_plus.append(v[:])

 #   print(f"----{Phi_0_plus}")
    #这里是所有奇正根
    Phi_1_plus = []

    for i in range(n):
        for j in range(m):
            v = zero_vector(QQ,n+m)
            v[i] = 1
            v[n+j] = 1
            Phi_1_plus.append(v[:])

    Coeff = zero_vector(QQ,n+m)
    for i in range(n):
        Coeff[i] = 1
    for i in range(m):
        Coeff[n+i] = -1
    
    lowest_module = Lowest_Module(n,m)
    Result = []
    Result_hash = set()


    lambda_beta = []
    lambda_init_class = Mu_With_Message(1,lambda_init)
    item = tuple(lambda_init)
    Result.append(lambda_init_class)
    Result_hash.add(item)

    lambda_beta.append(lambda_init_class)

    for beta in Phi_1_plus:
#        print(f"here---{lambda_init}*{beta}*{Coeff}={lambda_init*beta}*{Coeff} = {lambda_init * beta * Coeff}")
        if (lambda_init * vector(QQ,[ beta[i]*Coeff[i] for i in range(n+m)] ) ) == 0:
            item = tuple(lambda_init+beta)
            if item not in Result_hash:
                lambda_result = Mu_With_Message(2,lambda_init + beta)
                lambda_result.plus_roots.append(beta)
                lambda_beta.append(lambda_result)
                Result.append(lambda_result)
                Result_hash.add(item)

 #               print(f"1---{lambda_result.result}")

    for beta in Phi_1_plus:
        for gama in Phi_1_plus:
            flag_beta, Coe_beta =  is_nonnegative_integer_combination_sage(beta , lowest_module.basis_plus)
            flag_gama, Coe_gama =  is_nonnegative_integer_combination_sage(gama , lowest_module.basis_plus)
            
            if (lambda_init * vector(QQ,[ beta[i]*Coeff[i] for i in range(n+m)])) == 0 and ( (lambda_init + beta) * vector(QQ,[ gama[i]*Coeff[i] for i in range(n+m)]) ) == 0 and sum(Coe_beta) < sum(Coe_gama):

                item = tuple(lambda_init+beta+gama)
                if item not in Result_hash:
                    lambda_result = Mu_With_Message(3,lambda_init + beta+ gama)
                    lambda_result.plus_roots.append(beta)
                    lambda_result.plus_roots.append(gama)
                    lambda_beta.append(lambda_result)
                    Result.append(lambda_result)
                    Result_hash.add(item)
#                    print(f"2---{lambda_result.result}")
    

    for lambda_plus_beta_class in lambda_beta:

        lambda_plus_beta = lambda_plus_beta_class.result
        # 存储 (权重, 路径) 而不仅仅是权重
        lambda_before = []
        lambda_before.append((lambda_plus_beta[:], []))  # (权重, 路径)
        depth_1 = depth
#        if lambda_plus_beta_class.type == 1  :
#            print("测试开始")
#        else:
#            print("测试结束")
        while depth_1:
            depth_1 = depth_1 - 1
            lambda_after = []
            for lambda_con, path in lambda_before:
                for a in Phi_0_plus:
                    u = 2*(lambda_con * vector(QQ,[a[i]*Coeff[i] for i in range(n+m)]))/(a * vector(QQ,[ a[i]*Coeff[i] for i in range(n+m)]))
#                    print(f"{a}--{lambda_con}--{path}")
#                    print(vector(QQ,[a[i]*Coeff[i] for i in range(n+m)]))
#                    print(a * vector(QQ,[ a[i]*Coeff[i] for i in range(n+m)]))
#                    print(u)
                    if u < 0: 
                        lambda_new = lambda_con - u * a
                        item = tuple(lambda_new)
                        if item not in Result_hash:
                            # 记录新的路径：当前路径 + 使用的根向量
                            new_path = path + [tuple(a)]  # 存储a的副本

                            lambda_result = Mu_With_Message(lambda_plus_beta_class.type,lambda_new[:])
                            lambda_result.plus_roots = lambda_plus_beta_class.plus_roots
                            lambda_result.step_set = new_path
                            Result.append(lambda_result)
                            Result_hash.add(item)
                            lambda_after.append((lambda_new[:], new_path))

#                            print(f"------{lambda_result.result}")
            lambda_before.clear()
            lambda_before = lambda_after

    return Result_hash, Result





def show_steps(v,lambda_judge):
    for mu_class in lambda_judge:
        if mu_class.result == v:
            print(f"类型{mu_class.type}")
            if mu_class.type ==2 or mu_class.type ==3:
                print(f"路径{mu_class.step_set}")
                print(f"加权{mu_class.plus_roots}")
            else:
                print(f"路径{mu_class.step_set}")





def vectors_set_min(vector_set,vectors_set_min):
    vector_set_tem = vector_set[:]
    for v in vectors_set_min:
        try:
            vector_set_tem.remove(v)
        except ValueError:
            print(f"元素 {v} 不在列表中，无法移除")
            return None
    return vector_set_tem


def calc_typi_vec(typical_lambda_sp_plus_so,n,m):

    t_lambda_sp = typical_lambda_sp_plus_so[:n]
    t_lambda_so = typical_lambda_sp_plus_so[n:]
    print(f"typical的lambda是: {typical_lambda_sp_plus_so}")
    
    W_sp = Weyl_Group_Bn(n)
    W_so = Weyl_Group_Bn(m)
    
    w_sp, lambda_sp_next= calc_w_mu(W_sp,t_lambda_sp)
    w_so, lambda_so_next= calc_w_mu(W_so,-t_lambda_so)

    sum_sp_weyl,sum_so_weyl,sum_sp_plus_so = K_L_decompose(W_sp,w_sp,lambda_sp_next, W_so,w_so,lambda_so_next)  
    return sum_sp_weyl,sum_so_weyl,sum_sp_plus_so
