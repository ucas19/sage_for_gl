
from os import remove
from typing import Counter
#from functions import show_steps
#from functions import is_tensor_V_true
#from functions import is_tensor_V_true_not_show
from weyl_group_Bn import Weyl_Group_Bn
from sage_vector_store import SageVectorGroupStore
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.rings.rational_field import QQ
from lowest_module import Lowest_Module
from integer_com import is_nonnegative_integer_combination_simple
from sage_integer_com import is_nonnegative_integer_combination_sage
#from functions import *
from Mu_with_message import *
from functions_file import read_rational_vectors, read_vectors_from_file
from datetime import datetime
load("functions.sage")
def contains_with_counts(A, B):

    immutable_vecs_A = [vector(v, immutable=True) for v in A]
    immutable_vecs_B = [vector(v, immutable=True) for v in B]
    count_A = Counter(immutable_vecs_A)
    count_B = Counter(immutable_vecs_B)
    return all(count_A[x] >= count_B[x] for x in count_B)


def remove_matches_vec(typical_lambda_sp_plus_so,n,m):
    front = typical_lambda_sp_plus_so[:n] 
    back =  typical_lambda_sp_plus_so[n:]
    front_a, back_a, front_ctr, back_ctr, removed_count   =  remove_matches(front.list(),back.list())
    return removed_count

def add_store(lowest_weight,weight_set):

    print(f"即将储存特征标:{lowest_weight}")
    if store.exists(lowest_weight):
        items = store.get_group(lowest_weight)
        if not ( contains_with_counts(items,weight_set) and contains_with_counts(weight_set,items) ):
            current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            with open('a.txt', 'a', encoding='utf-8') as f:
                f.write("\n")
                f.write(f"-- {current_time} --警告,数据有问题")
                f.write("\n")
                f.write("\n")
            print(f"警告,数据有问题")
            print(f"警告,数据有问题")
            print(f"警告,数据有问题")
            print(f"警告,数据有问题")
            print(f"警告,数据有问题")

    try:
        store.add_group(lowest_weight, weight_set)
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        with open('a.txt', 'a', encoding='utf-8') as f:
            f.write("\n")
            f.write(f"-- {current_time} --写入特征标:{lowest_weight}")
            f.write(f"\n")
            f.write(f"{weight_set}")
            f.write(f"\n")
            f.write(f"\n")
    except ValueError as e:
        print(e)  # Key vector ... already exists



def pre_treatment_wheather_in_P(P_mu_tensor_V_after_Pr,n,m):
    print("")
    print(f"pre_treatment_wheather_in_P()")
    print("")
    user_sum = len(P_mu_tensor_V_after_Pr)# 将输入转换为有理数列表并创建向量
    not_consider_weight = []
    count = 1
    countss = 0
    in_consider_weight = []
    print(f"处理总数={user_sum}")
    zhengli(P_mu_tensor_V_after_Pr)

    immutable_vecs = [vector(v, immutable=True) for v in P_mu_tensor_V_after_Pr]
    P_mu_tensor_V_after_Pr_dic = Counter(immutable_vecs)
#    for v,k in P_mu_tensor_V_after_Pr_dic.items():
    for lambda_sp_plus_so, k in P_mu_tensor_V_after_Pr_dic.items():
        if store.exists(lambda_sp_plus_so):
            contains_lam_judge = []
            contains_lam_judge = store.get_group(lambda_sp_plus_so)
            print(" ")
            print(f"{count}: {lambda_sp_plus_so} 数量:{len(contains_lam_judge)}")
            print("--------------------")
            print(" ")
            count +=1
            if contains_with_counts(P_mu_tensor_V_after_Pr,contains_lam_judge):
                countss = countss+1
                in_consider_weight.append(lambda_sp_plus_so)
                print("在数据库找到了, 上面这个权可以要单独考虑")
                print("***********************")
            else:
                print("这个权不用考虑了。因为剩下的不能把这些直接计算的权包含进去")
            print("***********************")
            continue


        lambda_judge_hash, lambda_judge = judge_mu_in_P(lambda_sp_plus_so,n,m,20)

        if len(lambda_judge_hash) <= user_sum:
            print(" ")
            print(f"{count}: {lambda_sp_plus_so} 数量:{len(lambda_judge_hash)}")
            print("--------------------")
            print(" ")
            count +=1
            contains_lam_judge = []

            for item in lambda_judge:
                contains_lam_judge.append(item.result)

            if(contains_with_counts(P_mu_tensor_V_after_Pr,contains_lam_judge)):

                countss +=1
                count_count = 1
                for item in lambda_judge:
                    print(f"{count_count}: {item.result}")
                    count_count+=1
                in_consider_weight.append(lambda_sp_plus_so)
                print("上面这个权可以要单独考虑")
            else:
                print("这个权不用考虑了。因为剩下的不能把这些直接计算的权包含进去")
            print("***********************")
        else:
            print(" ")
            print(f"{count}: {lambda_sp_plus_so} 数量:{len(lambda_judge_hash)}")
            print("--------------------")
            print(" ")
            count+=1
            print("这个权不用考虑了。因为数量问题")
            not_consider_weight.append(lambda_sp_plus_so)

    print("以下因数量不需要考虑的向量是:")
    for i in range(len(not_consider_weight)):
        print(f"{i+1}: {not_consider_weight[i]}")
                
    print(" ")
    print(f"一共有{countss}个向量需要考虑")
    for i in range(len(in_consider_weight)):
        print(f"{i+1}: {in_consider_weight[i]}")

    return in_consider_weight 

 #           count = 1
 #           for item in lambda_judge:
 #               print(f"{count}: {item.result}")
 #               count+=1


def zhengli(sum_sp_plus_so):
    immutable_vecs = [vector(v, immutable=True) for v in sum_sp_plus_so]
    sum_sp_plus_so_count = Counter(immutable_vecs)
    count=1
    print("----------统计如下----------")
    for v,n in sum_sp_plus_so_count.items():
        print(f"{count}: {v} 数量{n}")
        count +=1
    print("----------整理结束----------")


def test_K_L(nn,mm,L_sp,L_so,flag=0):
    n = nn
    m = mm
    W_sp = Weyl_Group_Bn(n)
    W_so = Weyl_Group_Bn(m)
    
    lambda_sp = vector(QQ,L_sp)
    lambda_so = vector(QQ,L_so)
#    lambda_sp_plus_so = vector(QQ, L_sp_so)
#    print(f"lambda是: {lambda_sp_plus_so}")
    print(f"mu是: {vector(QQ,list(lambda_sp)+list(lambda_so))}")
    
    w_sp, lambda_sp_next= calc_w_mu(W_sp,lambda_sp)
    w_so, lambda_so_next= calc_w_mu(W_so,lambda_so)
    print(f"前下标weyl元是w_sp = {w_sp}")
    print(f"后下标weyl元是w_so = {w_so}")
    print(f"下标:{list(lambda_sp_next)+list(lambda_so_next)}")
    
    if flag == 0:
        sum_sp_weyl,sum_so_weyl,sum_sp_plus_so = K_L_decompose_no_kl(W_sp,w_sp,lambda_sp_next, W_so,w_so,lambda_so_next)  
        for i in range(len(sum_sp_plus_so)):
            print(f"{i+1} : {sum_sp_plus_so[i]}")

    else:
        sum_sp_weyl,sum_so_weyl,sum_sp_plus_so = K_L_decompose(W_sp,w_sp,lambda_sp_next, W_so,w_so,lambda_so_next)  
        immutable_vecs = [vector(v, immutable=True) for v in sum_sp_plus_so]
        sum_sp_plus_so_count = Counter(immutable_vecs)
        count=1
        print("----------统计如下----------")
        for v,n in sum_sp_plus_so_count.items():
            print(f"{count}: {v} 数量{n}")
            count +=1
        print(f"总数: {len(sum_sp_plus_so)}")
        print("+-+-+-+-+-+-+-+-+")
        print("展开:")
        print(sum_sp_plus_so)
    return sum_sp_plus_so
def show_kl_comps(P_mu_tensor_V_after_Pr,n,m):
    lowest_module = Lowest_Module(n,m)
    weight_set = P_mu_tensor_V_after_Pr[:]
    result = weight_set[:]
    result_low = []
    i=1
    while weight_set:
        lowest_weight = which_one_lowest(weight_set,lowest_module.basis_plus)
        lowest_weight_now = lowest_weight[0]
        lambda_sp = lowest_weight_now[:n]
        lambda_so = lowest_weight_now[n:]
        W_sp = Weyl_Group_Bn(n)
        W_so = Weyl_Group_Bn(m)
        w_sp, lambda_sp_next= calc_w_mu(W_sp,lambda_sp)
        w_so, lambda_so_next= calc_w_mu(W_so,lambda_so)
        sum_sp_weyl,sum_so_weyl,sum_sp_plus_so = K_L_decompose_no_kl(W_sp,w_sp,lambda_sp_next, W_so,w_so,lambda_so_next)  

        flag =1
        for v in sum_sp_plus_so:
            if v not in weight_set:
                flag = 0

        if flag:
            result_low.append(lowest_weight_now)
            print(f"第{i}个缩写: {lowest_weight_now}")
            i += 1
            for j in range(len(sum_sp_plus_so)):
                print(f"{j+1} : {sum_sp_plus_so[j]}")
                weight_set.remove(sum_sp_plus_so[j])
                result.remove(sum_sp_plus_so[j])
        else:
            weight_set.remove(lowest_weight_now)


    print("还剩下:")
    zhengli(result)







def again_calc(L_sp_so_next,P_after,which_mod,n,m):
    lambda_sp_plus_so = L_sp_so_next[:]
    sum_sp_plus_so = P_after[:]
    lowest_module = Lowest_Module(n,m)

    if which_mod ==1:
        P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr = P_tensor_V(lambda_sp_plus_so,sum_sp_plus_so,lowest_module.V,n,m)
    elif which_mod==2:
        P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr = P_tensor_V(lambda_sp_plus_so,sum_sp_plus_so,lowest_module.S2V,n,m)
    elif which_mod==3:
        P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr = P_tensor_V(lambda_sp_plus_so,sum_sp_plus_so,lowest_module.g,n,m)
    elif which_mod==4:
        P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr = P_tensor_V(lambda_sp_plus_so,sum_sp_plus_so,lowest_module.S3V,n,m)
    elif which_mod==5:
        P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr = P_tensor_V(lambda_sp_plus_so,sum_sp_plus_so,lowest_module.W3V,n,m)
    else:
        print("---------输入有误---------")
        

    flag = 0
    if is_tensor_V_true(lambda_sp_plus_so, P_mu_tensor_V_after_Pr, lowest_module.basis_plus):
        print(f"投射成立，lambda是极小权")
        flag = 1
    else:
        print(f"*********注意!*********")
        print(f"投射后不一定是最小权")
        return None

    lambda_judge_hash, lambda_judge = judge_mu_in_P(lambda_sp_plus_so,n,m,20)
    cant_judge = P_mu_tensor_V_after_Pr[:]
    if flag == 1:
#        lambda_judge_hash, lambda_judge = judge_mu_in_P(lambda_sp_plus_so,n,m,8)       
        check_counts = 0
        immutable_vecs = [vector(v, immutable=True) for v in P_mu_tensor_V_after_Pr]
        P_mu_tensor_V_after_Pr_dic = Counter(immutable_vecs)

        at_lambda_sp_plus_so_imm = vector(lambda_sp_plus_so, immutable=True)
        k_step = P_mu_tensor_V_after_Pr_dic[at_lambda_sp_plus_so_imm]

        for v,k in P_mu_tensor_V_after_Pr_dic.items():
            v_hash = tuple(v)
            if v_hash not in lambda_judge_hash:
                print(f"{v}不能判断是否在里面,数量: {k}")
            else:
                check_counts = check_counts+1
                for i in range(k_step):
                    cant_judge.remove(v)
#                print(f"{v}在里面")
#                show_steps(v,lambda_judge)
        if not check_counts == len(lambda_judge_hash):
            print("***********************")
            print("数据有问题, 本次计算作废!!!")
            print("***********************")
            flag = 0

        else:
            print("\n")
            print("-----计算kl折叠:------")
            show_kl_comps(P_mu_tensor_V_after_Pr,n,m)
#            else:
#                print(f"{v}在里面")
#                show_steps(v,lambda_judge)

    calc_sum = 1
    print("直接计算得到的在不在里面:")
    for result_ju in lambda_judge:
        print(f"{calc_sum}: {result_ju.result}")
        calc_sum +=1

    if flag ==1:

        if cant_judge:
            is_pre_treatment = input(f"是否需要处理本次结果")
            if is_pre_treatment == "yes":
                pre_treatment_wheather_in_P(cant_judge,n,m)

        is_store = input(f"是否储存本次特征标计算结果？")
        if is_store == "yes":
            add_store(lambda_sp_plus_so,P_mu_tensor_V_after_Pr)

    return P_mu_tensor_V_after_Pr


def again_calc_for_ten(L_sp_so_next,P_after,which_mod,n,m):
    print("")
    print(f"again_calc_for_ten({L_sp_so_next},which_mod = {which_mod})")
    print("")
    lambda_sp_plus_so = L_sp_so_next[:]
    sum_sp_plus_so = P_after[:]
    lowest_module = Lowest_Module(n,m)

    if which_mod ==1:
        P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr = P_tensor_V(lambda_sp_plus_so,sum_sp_plus_so,lowest_module.V,n,m)
    elif which_mod==2:
        P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr = P_tensor_V(lambda_sp_plus_so,sum_sp_plus_so,lowest_module.S2V,n,m)
    elif which_mod==3:
        P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr = P_tensor_V(lambda_sp_plus_so,sum_sp_plus_so,lowest_module.g,n,m)
    elif which_mod==4:
        P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr = P_tensor_V(lambda_sp_plus_so,sum_sp_plus_so,lowest_module.S3V,n,m)
    elif which_mod==5:
        P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr = P_tensor_V(lambda_sp_plus_so,sum_sp_plus_so,lowest_module.W3V,n,m)
    else:
        print("---------输入有误---------")
        

    flag = 0
    if is_tensor_V_true(lambda_sp_plus_so, P_mu_tensor_V_after_Pr, lowest_module.basis_plus):
        print(f"投射成立，lambda是极小权")
        flag = 1
    else:
        print(f"*********注意!*********")
        print(f"投射后不一定是最小权")
        return None

    return P_mu_tensor_V_after_Pr



def P_tensor_V_and_judge_wheather_minest(n,m,at_lambda_sp_plus_so,sum_sp_plus_so,which_mod):
    print("")
    print(f"P_tensor_V_and_judge_wheather_minest({at_lambda_sp_plus_so},which_mod = {which_mod})")
    print("")

    lowest_module = Lowest_Module(n,m)

    if which_mod ==1:
        P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr = P_tensor_V(at_lambda_sp_plus_so,sum_sp_plus_so,lowest_module.V,n,m)
    elif which_mod==2:
        P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr = P_tensor_V(at_lambda_sp_plus_so,sum_sp_plus_so,lowest_module.S2V,n,m)
    elif which_mod==3:
        P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr = P_tensor_V(at_lambda_sp_plus_so,sum_sp_plus_so,lowest_module.g,n,m)
    elif which_mod==4:
        P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr = P_tensor_V(at_lambda_sp_plus_so,sum_sp_plus_so,lowest_module.S3V,n,m)
    elif which_mod==5:
        P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr = P_tensor_V(at_lambda_sp_plus_so,sum_sp_plus_so,lowest_module.W3V,n,m)
    else:
        print("--------输入错误----------")
    flag = 0
    print(f"P_mu_tensor_V_befor_Pr: {len(P_mu_tensor_V_befor_Pr)}")
    print(f"P_mu_tensor_V_after_Pr: {len(P_mu_tensor_V_after_Pr)}")

    if is_tensor_V_true(at_lambda_sp_plus_so, P_mu_tensor_V_after_Pr, lowest_module.basis_plus):
        print(f"投射成立，lambda{at_lambda_sp_plus_so}是极小权")
        flag=1
    else:
        print(f"*********注意!*********")
        print(f"投射后不一定是最小权")
    return flag, P_mu_tensor_V_after_Pr


def judge_wheather_in_P(at_lambda_sp_plus_so, P_mu_tensor_V_after_Pr, lambda_judge_hash, lambda_judge, n, m):
    print("")
    print(f"judge_wheather_in_P({at_lambda_sp_plus_so})")
    print("")
    

    cant_judge = P_mu_tensor_V_after_Pr[:]
    check_counts = 0
    immutable_vecs = [vector(v, immutable=True) for v in P_mu_tensor_V_after_Pr]
    P_mu_tensor_V_after_Pr_dic = Counter(immutable_vecs)
    at_lambda_sp_plus_so_imm = vector(at_lambda_sp_plus_so, immutable=True)
    
    k_step = P_mu_tensor_V_after_Pr_dic[at_lambda_sp_plus_so_imm]
    print(f"结果{at_lambda_sp_plus_so}存在数量:{k_step}")


    flag = 1
    for v,k in P_mu_tensor_V_after_Pr_dic.items():
        v_hash = tuple(v)
        if v_hash not in lambda_judge_hash:
            print(f"{v}不能判断是否在里面,数量: {k}")
        else:
            check_counts = check_counts+1
            if k < k_step:
                print("***********************")
                print("数据有问题, 本次计算作废!!!")
                print("***********************")
                flag = 0
            for i in range(k_step):
                cant_judge.remove(v)
#            print(f"{v}在里面")
#            show_steps(v,lambda_judge)
    if not check_counts == len(lambda_judge_hash):
        print("***********************")
        print("数据有问题, 本次计算作废!!!")
        print("***********************")
        flag = 0

    else:
        print("\n")
        print("-----计算kl折叠:------")
        show_kl_comps(P_mu_tensor_V_after_Pr,n,m)

    return flag, cant_judge


def calc_direct_all(lambda_judge):
    calc_sum = 1
    print("直接计算得到的在不在里面:")
    for result_ju in lambda_judge:
        print(f"{calc_sum}: {result_ju.result}")
        calc_sum +=1


def deal_with_cant_judge(cant_judge, n, m):
    print("")
    print("deal_with_cant_judge")
    print("")

    need_deal_with_vec = pre_treatment_wheather_in_P(cant_judge,n,m)
    is_store = input(f"是否储存本次特征标计算结果？")
    if is_store == "yes":
        add_store(at_lambda_sp_plus_so,P_mu_tensor_V_after_Pr)


def circle_calc(all_vecs, need_deal_with_vec,which_mod,n,m,sheaf):
    print("")
    print(f"circle_calc(which_mod = {which_mod},sheaf = {sheaf})")
    print("")
    front_vecs = all_vecs[:]
#    at_calc_vec(nn,mm,typical_lambda_sp_plus_so,at_lambda_sp_plus_so,which_mod)

#    print(f"aaaa---{front}")

    select_blocks = P_tensor_V_show(front_vecs,which_mod,n,m)

    for deal_vec in need_deal_with_vec:
        if not store.exists(deal_vec):
            flag = auto_calc_long_time(store,deal_vec,n,m,sheaf) 
            if not flag:
                return flag

    flag = 1
#   print(f"xxxx--{front}")
    for deal_vec in need_deal_with_vec:
        behind_vecs = store.get_group(deal_vec)
#        print(f"xxxx--{front}")
#        print(f"yyyy--{behind}")
        P_weights = vectors_set_min(front_vecs,behind_vecs)
        if P_weights == None:
            print(f"P_weights中不能包含P_{deal_vec}")
            print("")
            return 1
        results_ten = []

        for key,groups in select_blocks.items():
#            print(f"===={n}===key={groups[0]}")
            
            front = groups[0][:n]
            back = groups[0][n:]
            front, back,front_ctr,back_ctr,removed_count = remove_matches(front.list(),back.list())
            #补充，这里可以选用递归计算，时间太久
            if not removed_count == 0:
                continue
                immutable_vecs = [vector(v, immutable=True) for v in groups]
                groups_dic = Counter(immutable_vecs)

                print(f"which_mod: {which_mod}")
                if which_mod ==1:
                    print("使用模:V")# 将输入转换为有理数列表并创建向量
                elif which_mod==2:
                    print("使用模:S2V")# 将输入转换为有理数列表并创建向量
                elif which_mod==3:
                    print("使用模:g")# 将输入转换为有理数列表并创建向量
                elif which_mod==4:
                    print("使用模:S3V")# 将输入转换为有理数列表并创建向量
                elif which_mod==5:
                    print("使用模:W3V")# 将输入转换为有理数列表并创建向量
                
                results_ten = []
                for vec, sums in groups_dic.items():                  
                    flag_2 = deal_with_atypi_ten(vec,P_weights,which_mod,n,m,sheaf)
                    if flag_2:
                        results_ten.append(vec)


            if removed_count == 0:
                immutable_vecs = [vector(v, immutable=True) for v in groups]
                groups_dic = Counter(immutable_vecs)
                print(f"which_mod: {which_mod}")
                if which_mod ==1:
                    print("使用模:V")# 将输入转换为有理数列表并创建向量
                elif which_mod==2:
                    print("使用模:S2V")# 将输入转换为有理数列表并创建向量
                elif which_mod==3:
                    print("使用模:g")# 将输入转换为有理数列表并创建向量
                elif which_mod==4:
                    print("使用模:S3V")# 将输入转换为有理数列表并创建向量
                elif which_mod==5:
                    print("使用模:W3V")# 将输入转换为有理数列表并创建向量
                
                results_ten = []
                for vec, sums in groups_dic.items():                  
                    flag_2 = deal_with_typi_ten(vec,P_weights,which_mod,n,m)
                    if flag_2:
                        results_ten.append(vec)


                print("--------------")
                print("打印结果")
                if results_ten==[]:
                    print("判断失效，都不行")
                    flag = 0
                else:
                    for i in range(len(results_ten)):
                        print(f"{i+1}: {results_ten[i]}")

    return flag







def test_a(nn,mm,typical_lambda_sp,typical_lambda_so,atypical_lambda_sp_plus_so,which_mod):

    n = nn
    m = mm
    W_sp = Weyl_Group_Bn(n)
    W_so = Weyl_Group_Bn(m)
    
    t_lambda_sp = vector(QQ,typical_lambda_sp)
    t_lambda_so = vector(QQ,typical_lambda_so)
    at_lambda_sp_plus_so = vector(QQ, atypical_lambda_sp_plus_so)
    print(f"atypical的lambda是: {at_lambda_sp_plus_so}")
    print(f"typical的lambda是: {vector(QQ,list(t_lambda_sp)+list(t_lambda_so))}")
    
    w_sp, lambda_sp_next= calc_w_mu(W_sp,t_lambda_sp)
    w_so, lambda_so_next= calc_w_mu(W_so,t_lambda_so)
    print(f"前下标weyl元是w_sp = {w_sp}")
    print(f"后下标weyl元是w_so = {w_so}")
    print(f"下标:{list(lambda_sp_next)+list(lambda_so_next)}")

    #    Dic_P_mu_tensor_V_befor_Pr = Counter(P_mu_tensor_V_befor_Pr)
#    Dic_P_mu_tensor_V_after_Pr = Counter(P_mu_tensor_V_after_Pr)

    sum_sp_weyl,sum_so_weyl,sum_sp_plus_so = K_L_decompose(W_sp,w_sp,lambda_sp_next, W_so,w_so,lambda_so_next)  

    lowest_module = Lowest_Module(n,m)
    print(f"sum_sp_plus_so: {len(sum_sp_plus_so)}")
    if which_mod ==1:
        P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr = P_tensor_V(at_lambda_sp_plus_so,sum_sp_plus_so,lowest_module.V,n,m)
    elif which_mod==2:
        P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr = P_tensor_V(at_lambda_sp_plus_so,sum_sp_plus_so,lowest_module.S2V,n,m)
    elif which_mod==3:
        P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr = P_tensor_V(at_lambda_sp_plus_so,sum_sp_plus_so,lowest_module.g,n,m)
    elif which_mod==4:
        P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr = P_tensor_V(at_lambda_sp_plus_so,sum_sp_plus_so,lowest_module.S3V,n,m)
    elif which_mod==5:
        P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr = P_tensor_V(at_lambda_sp_plus_so,sum_sp_plus_so,lowest_module.W3V,n,m)
    else:
        print("--------输入错误----------")
    flag = 0

    print(f"P_mu_tensor_V_befor_Pr: {len(P_mu_tensor_V_befor_Pr)}")
    print(f"P_mu_tensor_V_after_Pr: {len(P_mu_tensor_V_after_Pr)}")

    if is_tensor_V_true(at_lambda_sp_plus_so, P_mu_tensor_V_after_Pr, lowest_module.basis_plus):
        print(f"投射成立，lambda{at_lambda_sp_plus_so}是极小权")
        flag=1
    else:
        print(f"*********注意!*********")
        print(f"投射后不一定是最小权")

    lambda_judge_hash, lambda_judge = judge_mu_in_P(at_lambda_sp_plus_so,n,m,20)
    cant_judge = P_mu_tensor_V_after_Pr[:]
    if flag == 1:
#        lambda_judge_hash, lambda_judge = judge_mu_in_P(at_lambda_sp_plus_so,n,m,8)
        check_counts = 0
        immutable_vecs = [vector(v, immutable=True) for v in P_mu_tensor_V_after_Pr]
        P_mu_tensor_V_after_Pr_dic = Counter(immutable_vecs)
        for v,k in P_mu_tensor_V_after_Pr_dic.items():
            v_hash = tuple(v)
            if v_hash not in lambda_judge_hash:
                print(f"{v}不能判断是否在里面,数量: {k}")
            else:
                check_counts = check_counts+1
                cant_judge.remove(v)
#                print(f"{v}在里面")
#                show_steps(v,lambda_judge)
        if not check_counts == len(lambda_judge_hash):
            print("***********************")
            print("数据有问题, 本次计算作废!!!")
            print("***********************")
            flag = 0

        else:
            print("\n")
            print("-----计算kl折叠:------")
            show_kl_comps(P_mu_tensor_V_after_Pr,n,m)

    calc_sum = 1
    print("直接计算得到的在不在里面:")
    for result_ju in lambda_judge:
        print(f"{calc_sum}: {result_ju.result}")
        calc_sum +=1

    if flag == 1:
        if cant_judge:
            is_pre_treatment = input(f"是否需要处理本次结果")
            if is_pre_treatment == "yes":
                pre_treatment_wheather_in_P(cant_judge,n,m)
        is_store = input(f"是否储存本次特征标计算结果？")
        if is_store == "yes":
            add_store(at_lambda_sp_plus_so,P_mu_tensor_V_after_Pr)



def find_path_vector( lam, n, m , which_mod):

    lam_key = vector(QQ,[which_mod]+lam.list())
    if store_path.exists(lam_key):

        print("")
        print(f"正在进入find_path_vector({lam}),以前计算过,返回结果")
        print("")
        return store_path.get_group(lam_key)

    print("")
    print(f"正在进入find_path_vector({lam}),耗时比较长")
    print("")

    
    lowest_module = Lowest_Module(n,m)
    W_sp = Weyl_Group_Bn(n)
    W_so = Weyl_Group_Bn(m)

    lam_sp = lam[:n]
    lam_so = lam[n:]
    results_V = []
    results_S2V = []
    results_S3V = []
    results_W3V = []
    results_g = []
    if which_mod ==1:
        lowest_module_V = lowest_module.V
    if which_mod ==2:
        lowest_module_V = lowest_module.S2V
    if which_mod ==3:
        lowest_module_V = lowest_module.g
    if which_mod ==4:
        lowest_module_V = lowest_module.S3V
    if which_mod ==5:
        lowest_module_V = lowest_module.W3V


    anti_repeat_vectors_hash = []
    for w_sp in W_sp.W:
        for w_so in W_so.W:
            lam_sp_after = w_sp.to_matrix() * lam_sp
            lam_so_after = w_so.to_matrix() * lam_so
            sp_plus_so = vector(QQ, list(lam_sp_after)+list(lam_so_after))
            for v in lowest_module_V:
                at_lambda_sp_plus_so = v + sp_plus_so
                at_lambda_sp_plus_so_hash = tuple(at_lambda_sp_plus_so) 
                if at_lambda_sp_plus_so_hash in anti_repeat_vectors_hash:
                    continue
                anti_repeat_vectors_hash.append(at_lambda_sp_plus_so_hash)

    
                lambda_judge_hash, lambda_judge = judge_mu_in_P(at_lambda_sp_plus_so,n,m,8)
                P_vectors = []
                for result_ju in lambda_judge:
                    P_vectors.append(result_ju.result)
    
                P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr = P_tensor_V_not_show(lam, P_vectors,lowest_module_V,n,m)
                if is_tensor_V_true_not_show(lam,P_mu_tensor_V_after_Pr, lowest_module.basis_plus):
                    results_V.append(at_lambda_sp_plus_so)
    store_path.add_group(lam_key,results_V)
    return results_V
    """

    if which_mod ==2:
        for w_sp in W_sp.W:
            for w_so in W_so.W:
                lam_sp_after = w_sp.to_matrix() * lam_sp
                lam_so_after = w_so.to_matrix() * lam_so
                sp_plus_so = vector(QQ, list(lam_sp_after)+list(lam_so_after))
                for v in lowest_module.S2V:
                    at_lambda_sp_plus_so = v + sp_plus_so
 #                   print(f"test: {test}")
                    at_lambda_sp_plus_so_hash = tuple(at_lambda_sp_plus_so) 
                    if at_lambda_sp_plus_so_hash in anti_repeat_vectors_hash:
                        continue
                    anti_repeat_vectors_hash.append(at_lambda_sp_plus_so_hash)
    
                    lambda_judge_hash, lambda_judge = judge_mu_in_P(at_lambda_sp_plus_so,n,m,8)
                    P_vectors = []
                    for result_ju in lambda_judge:
                        P_vectors.append(result_ju.result)
    
                    P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr = P_tensor_V_not_show(lam, P_vectors,lowest_module.S2V,n,m)
                    if is_tensor_V_true_not_show(lam,P_mu_tensor_V_after_Pr, lowest_module.basis_plus):
                        results_S2V.append(at_lambda_sp_plus_so)
    if which_mod ==3:
        for w_sp in W_sp.W:
            for w_so in W_so.W:
                lam_sp_after = w_sp.to_matrix() * lam_sp
                lam_so_after = w_so.to_matrix() * lam_so
                sp_plus_so = vector(QQ, list(lam_sp_after)+list(lam_so_after))
                for v in lowest_module.g:
                    at_lambda_sp_plus_so = v + sp_plus_so
                    at_lambda_sp_plus_so_hash = tuple(at_lambda_sp_plus_so) 
                    if at_lambda_sp_plus_so_hash in anti_repeat_vectors_hash:
                        continue
                    anti_repeat_vectors_hash.append(at_lambda_sp_plus_so_hash)
    
                    lambda_judge_hash, lambda_judge = judge_mu_in_P(at_lambda_sp_plus_so,n,m,8)
                    P_vectors = []
                    for result_ju in lambda_judge:
                        P_vectors.append(result_ju.result)
    
                    P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr = P_tensor_V_not_show(lam, P_vectors,lowest_module.g,n,m)
                    if is_tensor_V_true_not_show(lam,P_mu_tensor_V_after_Pr, lowest_module.basis_plus):
                        results_g.append(at_lambda_sp_plus_so)

    if which_mod ==4:
        for w_sp in W_sp.W:
            for w_so in W_so.W:
                lam_sp_after = w_sp.to_matrix() * lam_sp
                lam_so_after = w_so.to_matrix() * lam_so
                sp_plus_so = vector(QQ, list(lam_sp_after)+list(lam_so_after))
                for v in lowest_module.S3V:
                    at_lambda_sp_plus_so = v + sp_plus_so
 #                   print(f"test: {test}")
                    at_lambda_sp_plus_so_hash = tuple(at_lambda_sp_plus_so) 
                    if at_lambda_sp_plus_so_hash in anti_repeat_vectors_hash:
                        continue
                    anti_repeat_vectors_hash.append(at_lambda_sp_plus_so_hash)
    
                    lambda_judge_hash, lambda_judge = judge_mu_in_P(at_lambda_sp_plus_so,n,m,8)
                    P_vectors = []
                    for result_ju in lambda_judge:
                        P_vectors.append(result_ju.result)
    
                    P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr = P_tensor_V_not_show(lam, P_vectors,lowest_module.S3V,n,m)
                    if is_tensor_V_true_not_show(lam,P_mu_tensor_V_after_Pr, lowest_module.basis_plus):
                        results_S3V.append(at_lambda_sp_plus_so)

    return results_V, results_S2V, results_g, results_S3V

    """

def deal_with_typi_ten(again_lam, P_weights, which_mod , n, m ):
        print("")
        print(f"deal_with_typi_ten( {again_lam}, P_weights,which_mod = {which_mod})")
        print("")

        P_mu_tensor_V_after_Pr = again_calc_for_ten(again_lam, P_weights, which_mod,n,m)

        if P_mu_tensor_V_after_Pr is None:
            print(f"不是极小权，结束循环")
            return 0
        print(f"循环要处理的总数: {len(P_mu_tensor_V_after_Pr)}")

 #       with open("test://front2.txt", "w", encoding="utf-8") as f:
 #           f.write("\n".join(P_mu_tensor_V_after_Pr))  # 每个元素一行

        weight_set = P_mu_tensor_V_after_Pr[:]
        lowest_module = Lowest_Module(n,m)
        flag = False
        while True:
            minest_tem = which_one_lowest( weight_set,lowest_module.basis_plus)
            for typical_lambda_sp_plus_so in minest_tem:
                print("------------------------")
                print(f"处理的极小权为: {typical_lambda_sp_plus_so} \n")
                print(f"剩下: {len(weight_set)}")
                typical_lambda_sp = typical_lambda_sp_plus_so[:n]
                typical_lambda_so = typical_lambda_sp_plus_so[-m:]
                
                results = test_K_L(n,m,typical_lambda_sp,typical_lambda_so,1)
                weight_set = vectors_set_min( weight_set,results ) 
                if weight_set is None:
                    break

            if weight_set is None:
                print(f"---- 恭喜你，结果是None ! ---- ")
                flag =True
                break

            if weight_set == []:
                print(f"---- 判断失效，请选择向量重新判断 ----")
                flag = False
                break
        return flag


def deal_with_atypi_ten(again_lam, P_weights, which_mod , n, m , sheaf ):
    print("")
    print(f"deal_with_atypi_ten( {again_lam}, P_weights,which_mod = {which_mod})")
    print("")

    P_mu_tensor_V_after_Pr = again_calc_for_ten(again_lam, P_weights, which_mod,n,m)

    if P_mu_tensor_V_after_Pr is None:
        print(f"不是极小权，结束循环")
        return 0
    print(f"循环要处理的总数: {len(P_mu_tensor_V_after_Pr)}")

 #   with open("test://front2.txt", "w", encoding="utf-8") as f:
 #       f.write("\n".join(P_mu_tensor_V_after_Pr))  # 每个元素一行

    weight_set = P_mu_tensor_V_after_Pr[:]
    lowest_module = Lowest_Module(n,m)
    flag = False
    while True:
        minest_tem = which_one_lowest( weight_set,lowest_module.basis_plus)
        for atypical_lambda_sp_plus_so in minest_tem:
            print("------------------------")
            print(f"处理的极小权为: {atypical_lambda_sp_plus_so} \n")
            print(f"剩下: {len(weight_set)}")

            if store.exists(atypical_lambda_sp_plus_so):
                results = store.get_group(atypical_lambda_sp_plus_so)
            else:
                if auto_calc_long_time(store,atypical_lambda_sp_plus_so,n,m,sheaf):
                    results = store.get_group(atypical_lambda_sp_plus_so)
                else:
                    return 0



            weight_set = vectors_set_min( weight_set,results ) 
            if weight_set is None:
                break

            if weight_set == []:
                break

        if weight_set is None:
            print(f"---- 恭喜你，结果是None ! ---- ")
            flag =True
            break

        if weight_set == []:
            print(f"---- 判断失效，请选择向量重新判断 ----")
            flag = False
            break
    return flag



def auto_calc_long_time(store, at_lambda_sp_plus_so,n,m, sheaf):
    print("")
    print(f"正在进入auto_calc_long_time({at_lambda_sp_plus_so},sheaf = {sheaf})")
    print("")

    if store.exists(at_lambda_sp_plus_so):
        print("已经存在")
        zhengli(store.get_group(at_lambda_sp_plus_so))
        return 1

    if not remove_matches_vec(at_lambda_sp_plus_so,n,m):
        print(f"{at_lambda_sp_plus_so}是typical的向量,结束！")
        return 0
    #补充，循环层数
    if sheaf == 3:
        return 0
    sheaf = sheaf+1

    #首先找到那些容易计算的向量
    wheather_end = 0
    """
    补充，这里可以拓展使用的总的投射模
    """
    for i in range(3):
        which_mod = i+1

        if wheather_end ==1:
            break

 #       alternative_vectors = read_vectors_from_file("test//"+"test2"+".txt")
        alternative_vectors = find_path_vector( at_lambda_sp_plus_so, n, m , which_mod)

        count = 1
        for v in alternative_vectors:
            print(f"{count}向量: {v}")
            count+=1
            

        lambda_judge_hash, lambda_judge = judge_mu_in_P(at_lambda_sp_plus_so,n,m,20)

        for typical_lambda_sp_plus_so in alternative_vectors:
            print("")
            print(f"使用{typical_lambda_sp_plus_so}计算{at_lambda_sp_plus_so}")
            print("")
            front = typical_lambda_sp_plus_so[:n] 
            back =  typical_lambda_sp_plus_so[n:]
            front_a, back_a, front_ctr, back_ctr, removed_count   =  remove_matches(front.list(),back.list())
            flag = 0
            P_weights = []
            if removed_count == 0:
                
                sum_sp_weyl,sum_so_weyl,P_weights = calc_typi_vec(typical_lambda_sp_plus_so,n,m)
            elif store.exists(typical_lambda_sp_plus_so):
                P_weights = store.get_group(typical_lambda_sp_plus_so)
            else:
                flag_4 = auto_calc_long_time(store,typical_lambda_sp_plus_so,n,m,sheaf)
                if flag_4:
                    P_weights = store.get_group(typical_lambda_sp_plus_so)
                else :
                    print(f"迭代次数过多:{typical_lambda_sp_plus_so}")
                    print("暂时不用这个向量计算了或者日后手动改写")
                    continue

                
            flag, P_mu_tensor_V_after_Pr = P_tensor_V_and_judge_wheather_minest(n,m,at_lambda_sp_plus_so,P_weights,which_mod)

            if flag:
                flag_2, cant_judge = judge_wheather_in_P(at_lambda_sp_plus_so, P_mu_tensor_V_after_Pr, lambda_judge_hash, lambda_judge, n, m)
                if flag_2 and not cant_judge:
                    print("")
                    print("已经计算出来了")
                    print(f"atypical向量:{at_lambda_sp_plus_so}")
                    print("")
                    zhengli(P_mu_tensor_V_after_Pr)
                    print("\n")
                    print("-----计算kl折叠:------")
                    show_kl_comps(P_mu_tensor_V_after_Pr,n,m)
                    add_store(at_lambda_sp_plus_so,P_mu_tensor_V_after_Pr)
                    wheather_end =1
                    break
                    
                if flag_2:
                    in_consider_weight = pre_treatment_wheather_in_P(cant_judge,n,m)
                    if not in_consider_weight:
                        print("")
                        print("已经计算出来了")
                        print(f"atypical向量:{at_lambda_sp_plus_so}")
                        print("")
                        zhengli(P_mu_tensor_V_after_Pr)
                        print("\n")
                        print("-----计算kl折叠:------")
                        show_kl_comps(P_mu_tensor_V_after_Pr,n,m)
                        add_store(at_lambda_sp_plus_so,P_mu_tensor_V_after_Pr)
                        wheather_end =1
                        break


                    flag_3 = 0
                    """
                    补充，这里可以拓展验证时使用的模
                    """
                    for i in range(3): 
                        which_mod_now = i+1
                        flag_3 = circle_calc(P_mu_tensor_V_after_Pr,in_consider_weight,which_mod_now,n,m,sheaf)
                        if flag_3:
                            break
                    if flag_3:
                        print("")
                        print("已经计算出来了")
                        print(f"atypical向量:{at_lambda_sp_plus_so}")
                        print("")
                        zhengli(P_mu_tensor_V_after_Pr)
                        print("\n")
                        print("-----计算kl折叠:------")
                        show_kl_comps(P_mu_tensor_V_after_Pr,n,m)
                        add_store(at_lambda_sp_plus_so,P_mu_tensor_V_after_Pr)
                        wheather_end =1
                        break
                        
    return wheather_end







if __name__ == "__main__":

    lambda_sp = vector(QQ,[3/2, 3/2])
    lambda_so = vector(QQ,[1/2])
    lambda_sp_plus_so = vector(QQ, [1/2,1/2,1/2])
    lambda_sp_plus_so_next = vector(QQ,[-1/2,1/2,1/2])

    
    n=2 
    m=2
    store = SageVectorGroupStore('my_vectors.db')
    store_path = SageVectorGroupStore('path_test_db.db')
    while True:
        print("*******开始计算***********************************")
        print("你想要做什么？")
        print("0,测试")
        print("1,直接计算")
        print("2,计算K_L分解")
        print("3,知道P,计算P_tensor_W并投射到块lambda")
        print("4,计算集合极小权")
        print("5,计算两个向量集合相减")
        print("6,直接计算typical权")
        print("7,计算P_tensor_W,但是不投射,找到极小权")
        print("8,整理权集")
        print("9,批量初步判断那些需要那些多余的atypical权集")
        print("10,慎用,一条龙计算,出发点是一个front文件, 一个behind文件, 一个lam文件")
        print("11,慎用,暴力寻解,计算速度极慢.计算可能能用来计算的权")
        print("12,文件中的向量转化成latex代码")
        print("13,特征标数据库操作")
        print("14,计算自动机,耗时极长,并且不一定有结果,适合离开工位时使用")

        while True:
            try:
                select_case = int(input("你的选择是: "))
                break  # 如果成功转换为整数，跳出循环
            except ValueError:
                print("输入无效，请输入一个整数。")
        if select_case==0:

            user_input = input("请输入文档的名字:")# 将输入转换为有理数列表并创建向量
            if user_input=='':
                print("返回操作")
                continue

            with open('test//'+user_input+'.txt', 'r') as f:
                lines = f.readlines()

            rat_list = [QQ(x.strip()) for x in lines[0].split(',')]
            key1 = vector(QQ,rat_list)          

            flag = auto_calc_long_time(store,key1,n,m,0)
            if flag:
                print("计算成功!")

           # lowest_module = Lowest_Module(n,m)
           # print(f"测试s3v:{len(lowest_module.S3V)}")
           # print(f"测试s3vv:{len(lowest_module.S3VV)}")
           # print(f"测试w3v:{len(lowest_module.W3V)}")
           # print(f"测试s2v:{len(lowest_module.S2V)}")
           # print(f"测试s2vv:{len(lowest_module.S2VV)}")
           # print(f"测试gV:{len(lowest_module.gV)}")
           # print(f"测试g:{len(lowest_module.g)}")
           # tem = vectors_set_min(lowest_module.S3VV,lowest_module.S3V)
           # print(tem)
           # tem = vectors_set_min(lowest_module.S2V,lowest_module.S2VV)
           # print(tem)
           # tem = vectors_set_min(lowest_module.g,lowest_module.gV)
           # print(tem)
           # test_kl(n)

        if select_case==14:
            print("读取要计算的特征标向量")
            user_input = input("请输入权集合set所在文档的名字:")# 将输入转换为有理数列表并创建向量
            weight_set = read_vectors_from_file("test//"+user_input+".txt")
            print("注意计算开始:")

            for v in weight_set:
                if flag == auto_calc_long_time(v,n,m):
                    print("计算成功!")
                    zhengli(store.get_group(v))


        if select_case==13:

            print("1, 查询特征标")
            print("2, 储存特征标,文件中的向量组转化储存起来备用")
            print("3, 删除特征标,特征标输入有误")
            print("4, 列出全部特征标keys")

            while True:
                try:
                    sub_select_case = int(input("你的选择是: "))
                    break  # 如果成功转换为整数，跳出循环
                except ValueError:
                    print("输入无效，请输入一个整数。")

            if sub_select_case == 4:
                store.list_keys()

            elif sub_select_case == 3:

                user_input = input("请输入文档的名字:")# 将输入转换为有理数列表并创建向量
                if user_input=='':
                    print("返回操作")
                    continue

                with open('test//'+user_input+'.txt', 'r') as f:
                    lines = f.readlines()

                rat_list = [QQ(x.strip()) for x in lines[0].split(',')]
                key1 = vector(QQ,rat_list)          
                print(f"即将删除: {rat_list}")
                store.remove_group(key1)
                print(f"删除成功: {rat_list}")


            elif sub_select_case == 2 :
                user_input = input("请输入权集合set所在文档的名字:")# 将输入转换为有理数列表并创建向量
                if user_input=='':
                    print("返回操作")
                    continue
                weight_set = read_vectors_from_file("test//"+user_input+".txt")
                lowest_module = Lowest_Module(n,m)
                lowest_weight = which_one_lowest(weight_set,lowest_module.basis_plus)
                if not len(lowest_weight) == 1:
                    print("特征标权集有误！")
                    continue
                else:
                    print(f"即将储存特征标:{lowest_weight[0]}")
                    try:
                        store.add_group(lowest_weight[0], weight_set)
                    except ValueError as e:
                        print(e)  # Key vector ... already exists
            elif sub_select_case == 1:

                user_input = input("请输入文档的名字:")# 将输入转换为有理数列表并创建向量
                with open('test//'+user_input+'.txt', 'r') as f:
                    lines = f.readlines()

                rat_list = [QQ(x.strip()) for x in lines[0].split(',')]
                key1 = vector(QQ,rat_list)          
                
                if not store.exists(key1):
                    print(f"还未储存向量:{key1}")
                    continue

                retrieved = store.get_group(key1)
                print(f"向量: {key1}")
                print(retrieved)  # [ (1, 0), (0, 1, 1/2) ]
                zhengli(retrieved)


        if select_case==12:
            user_input = input("请输入权集合set所在文档的名字:")# 将输入转换为有理数列表并创建向量
            weight_set = read_vectors_from_file("test//"+user_input+".txt")

            immutable_vecs = [vector(v, immutable=True) for v in weight_set]
            weight_set_count = Counter(immutable_vecs)

            if len(weight_set)==0:
                continue

            selects = input("你是否选择替换(回车表示不替换): ")

            if selects =='':
                for v,k in weight_set_count.items():
                    str_s = []
                    for i in range(len(v)):
                        if v[i]> 0:
                            str_s.append("\\fracc{"+str(abs(v[i])*2)+"}")
                        else:
                            str_s.append("-\\fracc{"+str(abs(v[i])*2)+"}")

                    results = ''
                    for i in range(n):
                        results = results + str_s[i] + ','
                    results = results[:-1] + '|' 
                    for i in range(m):
                        results = results + str_s[n+i] + ','
                    if k > 1:
                        print( str(k)+"M_{"+results[:-1] +"}")
                    else:
                        print( "M_{"+results[:-1] +"}")

            else: 
                num_tem_dict = {}
                for i in range(len(weight_set[0])):
                    num_tem = abs(weight_set[0][i])
                    if num_tem not in num_tem_dict:
                        char_c = input(f"{num_tem}代表: ")
                        if char_c == '':
                            num_tem_dict[num_tem] = str(num_tem) 
                            num_tem_dict[-num_tem] = '-'+str(num_tem) 
                        else:
                            char_c = ''.join(char_c.split())
                            if char_c[0] == '+':
                                char_c = char_c[1:]

                            num_tem_dict[num_tem] = char_c 

                            if len(char_c)==1:
                                num_tem_dict[-num_tem] = '-'+char_c
                            elif len(char_c)==2 and char_c[0] == '-':
                                num_tem_dict[-num_tem] =  char_c[1:]
                            elif len(char_c)==3 and char_c[1] == '+':
                                num_tem_dict[-num_tem] = '-'+char_c[0]+'-'+char_c[2]
                            elif len(char_c)==3 and char_c[1] == '-':
                                num_tem_dict[-num_tem] = '-'+char_c[0]+'+'+char_c[2]
                            elif len(char_c)==4 and char_c[0] == '-' and char_c[2] == '+':
                                num_tem_dict[-num_tem] = char_c[1]+'-'+char_c[3]
                            elif len(char_c)==4 and char_c[0] == '-' and char_c[2] == '-':
                                num_tem_dict[-num_tem] = char_c[1]+'+'+char_c[3]
                
                for v,k in weight_set_count.items():
                    str_s = []
                    for i in range(len(v)):
                            str_s.append(num_tem_dict[v[i]])

                    results = ''
                    for i in range(n):
                        results = results + str_s[i] + ','
                    results = results[:-1] + '|' 
                    for i in range(m):
                        results = results + str_s[n+i] + ','
                    if k > 1:
                        print(str(k)+"M_{"+results[:-1] +"}")
                    else:
                        print("M_{"+results[:-1] +"}")
                        

                            






        if select_case==1:
            user_input = input("请输入文档的名字:")# 将输入转换为有理数列表并创建向量
            with open('test//'+user_input+'.txt', 'r') as f:
                lines = f.readlines()

            rat_list = [QQ(x.strip()) for x in lines[0].split(',')]
            atypical_lambda_sp_plus_so = vector(QQ,rat_list)          
            print(f"要计算的atipycal有理数向量lambda(用逗号,分隔):{atypical_lambda_sp_plus_so}")# 将输入转换为有理数列表并创建向量


            rat_list = [QQ(x.strip()) for x in lines[1].split(',')]
            typical_lambda_sp_plus_so = vector(QQ, rat_list)          
            
            print(f"要计算的tipycal有理数向量lambda(用逗号,分隔):{typical_lambda_sp_plus_so}")# 将输入转换为有理数列表并创建向量
            front = typical_lambda_sp_plus_so[:n] 
            back =  typical_lambda_sp_plus_so[n:]
            front_a, back_a, front_ctr, back_ctr, removed_count   =  remove_matches(front.list(),back.list())
#            print(f"---------{removed_count}")
            if not removed_count == 0:
                print(f"请检查, 并非typical向量:{typical_lambda_sp_plus_so}")
                continue

            which_mod = int(lines[2])
            if which_mod ==1:
                print("使用模:V")# 将输入转换为有理数列表并创建向量
            elif which_mod==2:
                print("使用模:S2V")# 将输入转换为有理数列表并创建向量
            elif which_mod==3:
                print("使用模:g")# 将输入转换为有理数列表并创建向量
            elif which_mod==4:
                print("使用模:S3V")# 将输入转换为有理数列表并创建向量
            elif which_mod==5:
                print("使用模:W3V")# 将输入转换为有理数列表并创建向量
            typical_lambda_sp = typical_lambda_sp_plus_so[:n]
            typical_lambda_so = typical_lambda_sp_plus_so[-m:]
            test_a(n,m,typical_lambda_sp,typical_lambda_so, atypical_lambda_sp_plus_so,which_mod)
        elif select_case==2:

            user_input = input("请输入文档的名字:")# 将输入转换为有理数列表并创建向量
            with open('test//'+user_input+'.txt', 'r') as f:
                lines = f.readlines()
            for line in lines:
                rat_list = [QQ(x.strip()) for x in line.split(',')]
                typical_lambda_sp_plus_so = vector(QQ,rat_list)          
                print(f"要计算的tipycal有理数向量lambda(用逗号,分隔):{typical_lambda_sp_plus_so}")# 将输入转换为有理数列表并创建向量
                typical_lambda_sp = typical_lambda_sp_plus_so[:n]
                typical_lambda_so = typical_lambda_sp_plus_so[-m:]
                test_K_L(n,m,typical_lambda_sp,typical_lambda_so)
        elif select_case==3:
                
            print("1.文件取用")
            print("2.数据库取用")
            while True:
                try:
                    sub_select_case = int(input("你的选择是: "))
                    break  # 如果成功转换为整数，跳出循环
                except ValueError:
                    print("输入无效，请输入一个整数。")
            if sub_select_case == 2:
                user_input = input("请输入文档的名字:")# 将输入转换为有理数列表并创建向量
                with open('test//'+user_input+'.txt', 'r') as f:
                    lines = f.readlines()

                rat_list = [QQ(x.strip()) for x in lines[0].split(',')]
                atypical_lambda_sp_plus_so = vector(QQ,rat_list)          
                print(f"要计算的atipycal有理数向量lambda(用逗号,分隔):{atypical_lambda_sp_plus_so}")# 将输入转换为有理数列表并创建向量


                rat_list = [QQ(x.strip()) for x in lines[1].split(',')]
                typical_lambda_sp_plus_so = vector(QQ, rat_list)          
                print(f"使用的数据库中有理数向量lambda(用逗号,分隔):{typical_lambda_sp_plus_so}")# 将输入转换为有理数列表并创建向量
                

                if not store.exists(typical_lambda_sp_plus_so) :
                    print("数据库中不存在这个特征标")
                    continue

                P_mu_tensor_V_after_Pr = store.get_group(typical_lambda_sp_plus_so)

                which_mod = int(lines[2])
                if which_mod ==1:
                    print("使用模:V")# 将输入转换为有理数列表并创建向量
                elif which_mod==2:
                    print("使用模:S2V")# 将输入转换为有理数列表并创建向量
                elif which_mod==3:
                    print("使用模:g")# 将输入转换为有理数列表并创建向量
                elif which_mod==4:
                    print("使用模:S3V")# 将输入转换为有理数列表并创建向量
                elif which_mod==5:
                    print("使用模:W3V")# 将输入转换为有理数列表并创建向量

                P_mu_tensor_V_after_Pr = again_calc(atypical_lambda_sp_plus_so,P_mu_tensor_V_after_Pr,which_mod,n,m)
            elif sub_select_case ==1:

                user_input = input("请输入P所在文档的名字:")# 将输入转换为有理数列表并创建向量
                    
                P_mu_tensor_V_after_Pr = read_vectors_from_file("test//"+user_input+".txt")
                
                user_input = input("请输入要投射的块lambda所在文档名字:")# 将输入转换为有理数列表并创建向量
                with open('test//'+user_input+'.txt', 'r') as f:
                    lines = f.readlines()
                rat_list = [QQ(x.strip()) for x in lines[0].split(',')]
                again_lam = vector(QQ,rat_list)          
                print(f"要投射的块是{again_lam}")

                which_mod = int(lines[1])
                if which_mod ==1:
                    print("使用模:V")# 将输入转换为有理数列表并创建向量
                elif which_mod==2:
                    print("使用模:S2V")# 将输入转换为有理数列表并创建向量
                elif which_mod==3:
                    print("使用模:g")# 将输入转换为有理数列表并创建向量
                elif which_mod==4:
                    print("使用模:S3V")# 将输入转换为有理数列表并创建向量
                elif which_mod==5:
                    print("使用模:W3V")# 将输入转换为有理数列表并创建向量
                P_mu_tensor_V_after_Pr = again_calc(again_lam,P_mu_tensor_V_after_Pr,which_mod,n,m)
        elif select_case==4:
            user_input = input("请输入权集合set所在文档的名字:")# 将输入转换为有理数列表并创建向量
            weight_set = read_vectors_from_file("test//"+user_input+".txt")
            lowest_module = Lowest_Module(n,m)
            print(f"1:完整权集大小排序")
            print(f"2:仅计算最小权")
            user_input = input("我的选择:")# 将输入转换为有理数列表并创建向量
            if user_input == '' or user_input =='2':
                lowest_weight = which_one_lowest(weight_set,lowest_module.basis_plus)
                print(f"极小权有{cases}: ")
                for v in lowest_weight:
                    print(f" {v} ")
                print("---------------------")
                continue
            cases = 1

            while weight_set:
                lowest_weight = which_one_lowest(weight_set,lowest_module.basis_plus)
                print(f"极小权有{cases}: ")
                for v in lowest_weight:
                    print(f" {v} ")
                print("---------------------")
                for v in lowest_weight:
                    weight_set.remove(v)

                cases = cases+1

        elif select_case==5:
            user_input_1 = input("请输入第一个权集合set所在文档的名字:")# 将输入转换为有理数列表并创建向量
            weight_set = read_vectors_from_file("test//"+user_input_1+".txt")
            user_input_2 = input("请输入第一个权集合set所在文档的名字:")# 将输入转换为有理数列表并创建向量
            weight_set_2 = read_vectors_from_file("test//"+user_input_2+".txt")
            weight_set = vectors_set_min(weight_set,weight_set_2)
            print("结果如下")
            for v in weight_set:
                print(v)
        elif select_case==6:
            user_input = input("请输入文档的名字:")# 将输入转换为有理数列表并创建向量
            with open('test//'+user_input+'.txt', 'r') as f:
                lines = f.readlines()
            for line in lines:
                rat_list = [QQ(x.strip()) for x in line.split(',')]
                typical_lambda_sp_plus_so = vector(QQ,rat_list)          
                print(f"直接计算typical情况，要计算的tipycal有理数向量lambda(用逗号,分隔):{typical_lambda_sp_plus_so}")# 将输入转换为有理数列表并创建向量
                typical_lambda_sp = typical_lambda_sp_plus_so[:n]
                typical_lambda_so = typical_lambda_sp_plus_so[-m:]
                test_K_L(n,m,typical_lambda_sp,typical_lambda_so,1)
        
        elif select_case==7:
                
            user_input = input("请输入P所在文档的名字:")# 将输入转换为有理数列表并创建向量
                
            P_mu_tensor_V_after_Pr = read_vectors_from_file("test//"+user_input+".txt")
            
            user_input_W = input("请输入模1，2，3，4，5:  ")# 将输入转换为有理数列表并创建向量

            which_mod = int(user_input_W)
            if which_mod ==1:
                print("使用模:V")# 将输入转换为有理数列表并创建向量
                P_tensor_V_show(P_mu_tensor_V_after_Pr,1,n,m)
            elif which_mod==2:
                print("使用模:S2V")# 将输入转换为有理数列表并创建向量
                P_tensor_V_show(P_mu_tensor_V_after_Pr,2,n,m)
            elif which_mod==3:
                print("使用模:g")# 将输入转换为有理数列表并创建向量
                P_tensor_V_show(P_mu_tensor_V_after_Pr,3,n,m)
            elif which_mod==4:
                print("使用模:S3V")# 将输入转换为有理数列表并创建向量
                P_tensor_V_show(P_mu_tensor_V_after_Pr,4,n,m)
            elif which_mod==5:
                print("使用模:W3V")# 将输入转换为有理数列表并创建向量
                P_tensor_V_show(P_mu_tensor_V_after_Pr,5,n,m)

        elif select_case==8:
            user_input = input("请输入权集所在文档的名字:")# 将输入转换为有理数列表并创建向量
                
            P_mu_tensor_V_after_Pr = read_vectors_from_file("test//"+user_input+".txt")
            zhengli(P_mu_tensor_V_after_Pr)
            show_kl_comps(P_mu_tensor_V_after_Pr,n,m)

        elif select_case==9:
                
            user_input = input("请输入P所在文档的名字:")# 将输入转换为有理数列表并创建向量
            
            results = []
            P_mu_tensor_V_after_Pr = read_vectors_from_file("test//"+user_input+".txt")
            user_sum = len(P_mu_tensor_V_after_Pr)# 将输入转换为有理数列表并创建向量
            
            not_consider_weight = []
            count = 1
            countss = 0
            in_consider_weight = []
            for lambda_sp_plus_so in P_mu_tensor_V_after_Pr:
                lambda_judge_hash, lambda_judge = judge_mu_in_P(lambda_sp_plus_so,n,m,20)

                if len(lambda_judge_hash) <= user_sum:
                    print(" ")
                    print(f"{count}: {lambda_sp_plus_so} 数量:{len(lambda_judge_hash)}")
                    print("--------------------")
                    print(" ")
                    count +=1
                    contains_lam_judge = []

                    for item in lambda_judge:
                        contains_lam_judge.append(item.result)

                    if(contains_with_counts(P_mu_tensor_V_after_Pr,contains_lam_judge)):

                        countss +=1
                        count_count = 1
                        for item in lambda_judge:
                            print(f"{count_count}: {item.result}")
                            count_count+=1
                        in_consider_weight.append(lambda_sp_plus_so)
                        print("上面这个权可以要单独考虑")
                    else:
                        print("这个权不用考虑了。因为剩下的不能把这些直接计算的权包含进去")
                    print("***********************")
                else:
                    print(" ")
                    print(f"{count}: {lambda_sp_plus_so} 数量:{len(lambda_judge_hash)}")
                    print("--------------------")
                    print(" ")
                    count+=1
                    print("这个权不用考虑了。因为数量问题")
                    not_consider_weight.append(lambda_sp_plus_so)

            print("以下因数量不需要考虑的向量是:")
            for i in range(len(not_consider_weight)):
                print(f"{i+1}: {not_consider_weight[i]}")
                        
            print(" ")
            print(f"一共有{countss}个向量需要考虑")
            for i in range(len(in_consider_weight)):
                print(f"{i+1}: {in_consider_weight[i]}")

 #               count = 1
 #               for item in lambda_judge:
 #                   print(f"{count}: {item.result}")
 #                   count+=1
        elif select_case == 10:


            user_input_front = input("请输入front(全部向量)所在文档的名字:")# 将输入转换为有理数列表并创建向量
            user_input_behind = input("请输入behind(要减去的向量)所在文档的名字:")# 将输入转换为有理数列表并创建向量
            user_input_lam = input("请输入lam(要投射的向量)所在文档的名字:")# 将输入转换为有理数列表并创建向量
            front = read_vectors_from_file("test//"+user_input_front+".txt")
            behind = read_vectors_from_file("test//"+user_input_behind+".txt")
            lam_s = read_vectors_from_file("test//"+user_input_lam+".txt")
            P_weights = vectors_set_min(front,behind)

#            with open('test//'+user_input_lam +'.txt', 'r') as f:
#                lines = f.readlines()
#            rat_list = [QQ(x.strip()) for x in lines[0].split(',')]
#            again_lam = vector(QQ,rat_list)          

            while True:
                try:
                    which_mod = int(input("你的选择是: "))
                    break  # 如果成功转换为整数，跳出循环
                except ValueError:
                    print("输入无效，请输入一个整数。")
            print(f"which_mod: {which_mod}")

#            which_mod = int(lines[1])
            if which_mod ==1:
                print("使用模:V")# 将输入转换为有理数列表并创建向量
            elif which_mod==2:
                print("使用模:S2V")# 将输入转换为有理数列表并创建向量
            elif which_mod==3:
                print("使用模:g")# 将输入转换为有理数列表并创建向量
            elif which_mod==4:
                print("使用模:S3V")# 将输入转换为有理数列表并创建向量
            elif which_mod==5:
                print("使用模:W3V")# 将输入转换为有理数列表并创建向量

            results_ten = []
            lam_s_imm = [vector(v, immutable=True) for v in lam_s]
            lam_s_imm_dict = Counter(lam_s_imm)
            for lam, kk in lam_s_imm_dict.items():
                if not remove_matches_vec(lam,n,m)==0:
                    continue

                flag = deal_with_typi_ten(lam,P_weights,which_mod,n,m)
                if flag:
                    results_ten.append(lam)


            print("--------------")
            print("打印结果")
            if results_ten==[]:
                print("判断失效，都不行")
            else:
                for i in range(len(results_ten)):
                    print(f"{i+1}: {results_ten[i]}")

#            P_mu_tensor_V_after_Pr = again_calc(again_lam, P_weights, which_mod,n,m)
#            if P_mu_tensor_V_after_Pr is None:
#                print(f"不是极小权，结束循环")
#                continue
#            print(f"循环要处理的总数: {len(P_mu_tensor_V_after_Pr)}")
#
# #           with open("test://front2.txt", "w", encoding="utf-8") as f:
# #               f.write("\n".join(P_mu_tensor_V_after_Pr))  # 每个元素一行
#
#            weight_set = P_mu_tensor_V_after_Pr[:]
#            lowest_module = Lowest_Module(n,m)
#            while True:
#                minest_tem = which_one_lowest( weight_set,lowest_module.basis_plus)
#                for typical_lambda_sp_plus_so in minest_tem:
#                    print("------------------------")
#                    print(f"处理的极小权为: {typical_lambda_sp_plus_so} \n")
#                    print(f"剩下: {len(weight_set)}")
#                    typical_lambda_sp = typical_lambda_sp_plus_so[:n]
#                    typical_lambda_so = typical_lambda_sp_plus_so[-m:]
#                    
#                    results = test_K_L(n,m,typical_lambda_sp,typical_lambda_so,1)
#                    weight_set = vectors_set_min( weight_set,results ) 
#                    if weight_set is None:
#                        break
#
#                if weight_set is None:
#                    print(f"---- 恭喜你，结果是None ! ---- ")
#                    break
#
#                if weight_set == []:
#                    print(f"---- 判断失效，请选择向量重新判断 ----")
#                    break

        elif select_case == 11:

            user_input = input("请输入文档的名字:")# 将输入转换为有理数列表并创建向量
            with open('test//'+user_input+'.txt', 'r') as f:
                lines = f.readlines()

            rat_list = [QQ(x.strip()) for x in lines[0].split(',')]
            atypical_lambda_sp_plus_so = vector(QQ,rat_list)          

            print(f"要计算的权是{atypical_lambda_sp_plus_so}")
            while True:
                try:
                    which_mod = int(input("你的选择是: "))
                    break  # 如果成功转换为整数，跳出循环
                except ValueError:
                    print("输入无效，请输入一个整数。")
            print(f"which_mod: {which_mod}")
            results_V = find_path_vector( atypical_lambda_sp_plus_so, n, m,which_mod )
            #results_V, results_S2V, results_g,results_S3V ,results_W3V = find_path_vector( atypical_lambda_sp_plus_so, n, m,which_mod )
            print(f"结果如下 V:{which_mod}")
            results = results_V

            """
            if which_mod ==1:
                results = results_V
            elif which_mod ==2:
                results = results_S2V
            elif which_mod ==3:
                results = results_g
            elif which_mod ==4:
                results = results_S3V
            elif which_mod ==4:
                results = results_W3V
            """

            count = 1
            for v in results:
                print(f"{count}向量: {v}")
                count+=1
            

        else:

            continue
    
    test_K_L(2,1,vector(QQ,[-3/2, -1/2]),vector(QQ,[3/2]))
    test_K_L(2,1,vector(QQ,[-1/2, -1/2]),vector(QQ,[1/2]))
    """
    again_lam = vector(QQ,[-3/2,-1/2,3/2])

    print(f"现在计算{again_lam}------------------")
    P_mu_tensor_V_after_Pr = read_vectors_from_file("test//test.txt")
    P_mu_tensor_V_after_Pr= again_calc(again_lam,P_mu_tensor_V_after_Pr,2,1)

    P_mu_tensor_V_after_Pr = test_a(2,1,lambda_sp,lambda_so,lambda_sp_plus_so,lambda_sp_plus_so_next)

    lambda_judge_hash, lambda_judge = judge_mu_in_P(lambda_sp_plus_so,2,1,5)
    for a in lambda_judge:
        print(a.result)
    
    again_lam = vector(QQ,[-3/2,1/2,1/2])
    print(f"现在计算{again_lam}------------------")
    P_mu_tensor_V_after_Pr= again_calc(again_lam,P_mu_tensor_V_after_Pr,2,1)
    
    again_lam = vector(QQ,[-3/2,-1/2,3/2])
    print(f"现在计算{again_lam}------------------")
    P_mu_tensor_V_after_Pr = read_vectors_from_file("test//n3_n1_1.txt")
    P_mu_tensor_V_after_Pr= again_calc(again_lam,P_mu_tensor_V_after_Pr,2,1)

    lowest_module = Lowest_Module(2,1)
    lowest_weight_set = read_vectors_from_file("test//test2.txt")
    lowest_weight = which_one_lowest(P_mu_tensor_V_after_Pr,lowest_module.basis_plus)
    print(f"极小权有{lowest_weight}")
    """

    lowest_module = Lowest_Module(2,1)
    weight_set = read_vectors_from_file("test//Pr_n3_n1_3.txt")
    weight_set_min = read_vectors_from_file("test//n3_n1_3.txt")
    weight_set = vectors_set_min(weight_set,weight_set_min)
    lowest_weight = which_one_lowest(weight_set,lowest_module.basis_plus)
    print(f"极小权有{lowest_weight}")
    
    weight_set_min = read_vectors_from_file("test//n3_1_n3.txt")
    weight_set = vectors_set_min(weight_set,weight_set_min)
    lowest_weight = which_one_lowest(weight_set,lowest_module.basis_plus)
    print(f"极小权有{lowest_weight}在下面")
    print(weight_set)

    #    test_K_L(2,1,vector(QQ,[-3/2, 1/2]),vector(QQ,[1/2]))
    #test_K_L(2,1,vector(QQ,[-3/2, 3/2]),vector(QQ,[3/2]))
    #test_K_L(2,1,vector(QQ,[-1/2, 3/2]),vector(QQ,[1/2]))


    #    lambda_judge_hash, lambda_judge = judge_mu_in_P(lambda_sp_plus_so,n,m,8)
    #for l_ju in lambda_judge:
    #    print(l_ju.result)

    again_lam = lowest_weight[0]
    print(f"现在计算{again_lam}------------------")
    lambda_judge_hash, lambda_judge = judge_mu_in_P(again_lam,2,1,20)
    for v in weight_set:
        v_hash = tuple(v)
        if v_hash not in lambda_judge_hash:
            print(f"{v}不能判断是否在里面")
        else:
            print(f"{v}在里面")
            show_steps(v,lambda_judge)

    weight_set = read_vectors_from_file("test//test3.txt")
    lowest_weight = which_one_lowest(weight_set,lowest_module.basis_plus)
    print(lowest_weight)
