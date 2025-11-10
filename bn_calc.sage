
from typing import Counter
#from functions import show_steps
from weyl_group_Bn import Weyl_Group_Bn
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.rings.rational_field import QQ
from lowest_module import Lowest_Module
from integer_com import is_nonnegative_integer_combination_simple
from sage_integer_com import is_nonnegative_integer_combination_sage
#from functions import *
from Mu_with_message import *
from functions_file import read_rational_vectors, read_vectors_from_file
load("functions.sage")
def contains_with_counts(A, B):

    immutable_vecs_A = [vector(v, immutable=True) for v in A]
    immutable_vecs_B = [vector(v, immutable=True) for v in B]
    count_A = Counter(immutable_vecs_A)
    count_B = Counter(immutable_vecs_B)
    return all(count_A[x] >= count_B[x] for x in count_B)

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
    result = []
    result_low = []
    flag =1
    i=1
    while weight_set and flag:
        result = weight_set[:]
        lowest_weight = which_one_lowest(weight_set,lowest_module.basis_plus)
        lowest_weight_now = lowest_weight[0]
        lambda_sp = lowest_weight_now[:n]
        lambda_so = lowest_weight_now[n:]
        W_sp = Weyl_Group_Bn(n)
        W_so = Weyl_Group_Bn(m)
        w_sp, lambda_sp_next= calc_w_mu(W_sp,lambda_sp)
        w_so, lambda_so_next= calc_w_mu(W_so,lambda_so)
        sum_sp_weyl,sum_so_weyl,sum_sp_plus_so = K_L_decompose_no_kl(W_sp,w_sp,lambda_sp_next, W_so,w_so,lambda_so_next)  
        for v in sum_sp_plus_so:
            if v in weight_set:
                weight_set.remove(v)
            else:
                flag = 0
        if flag:
            result_low.append(lowest_weight_now)
            print(f"第{i}个缩写: {lowest_weight_now}")
            i += 1
            for j in range(len(sum_sp_plus_so)):
                print(f"{j+1} : {sum_sp_plus_so[j]}")










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
    if flag == 1:
#        lambda_judge_hash, lambda_judge = judge_mu_in_P(lambda_sp_plus_so,n,m,8)
        immutable_vecs = [vector(v, immutable=True) for v in P_mu_tensor_V_after_Pr]
        P_mu_tensor_V_after_Pr_dic = Counter(immutable_vecs)
        for v,k in P_mu_tensor_V_after_Pr_dic.items():
            v_hash = tuple(v)
            if v_hash not in lambda_judge_hash:
                print(f"{v}不能判断是否在里面")
#            else:
#                print(f"{v}在里面")
#                show_steps(v,lambda_judge)
        print("\n")
        print("-----计算kl折叠:------")
        show_kl_comps(P_mu_tensor_V_after_Pr,n,m)

    calc_sum = 1
    print("直接计算得到的在不在里面:")
    for result_ju in lambda_judge:
        print(f"{calc_sum}: {result_ju.result}")
        calc_sum +=1
    return P_mu_tensor_V_after_Pr


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
    #    print(sum_sp_plus_so)
    if which_mod ==1:
        P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr = P_tensor_V(at_lambda_sp_plus_so,sum_sp_plus_so,lowest_module.V,n,m)
    elif which_mod==2:
        P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr = P_tensor_V(at_lambda_sp_plus_so,sum_sp_plus_so,lowest_module.S2V,n,m)
    elif which_mod==3:
        P_mu_tensor_V_befor_Pr, P_mu_tensor_V_after_Pr = P_tensor_V(at_lambda_sp_plus_so,sum_sp_plus_so,lowest_module.g,n,m)
    else:
        print("--------输入错误----------")
    flag = 0
    if is_tensor_V_true(at_lambda_sp_plus_so, P_mu_tensor_V_after_Pr, lowest_module.basis_plus):
        print(f"投射成立，lambda{at_lambda_sp_plus_so}是极小权")
        flag=1
    else:
        print(f"*********注意!*********")
        print(f"投射后不一定是最小权")

    lambda_judge_hash, lambda_judge = judge_mu_in_P(at_lambda_sp_plus_so,n,m,20)
    if flag == 1:
#        lambda_judge_hash, lambda_judge = judge_mu_in_P(at_lambda_sp_plus_so,n,m,8)

        immutable_vecs = [vector(v, immutable=True) for v in P_mu_tensor_V_after_Pr]
        P_mu_tensor_V_after_Pr_dic = Counter(immutable_vecs)
        for v,k in P_mu_tensor_V_after_Pr_dic.items():
            v_hash = tuple(v)
            if v_hash not in lambda_judge_hash:
                print(f"{v}不能判断是否在里面")
#            else:
#                print(f"{v}在里面")
#                show_steps(v,lambda_judge)
        print("\n")
        print("-----计算kl折叠:------")
        show_kl_comps(P_mu_tensor_V_after_Pr,n,m)

    calc_sum = 1
    print("直接计算得到的在不在里面:")
    for result_ju in lambda_judge:
        print(f"{calc_sum}: {result_ju.result}")
        calc_sum +=1




def zhengli(sum_sp_plus_so):
    immutable_vecs = [vector(v, immutable=True) for v in sum_sp_plus_so]
    sum_sp_plus_so_count = Counter(immutable_vecs)
    count=1
    print("----------统计如下----------")
    for v,n in sum_sp_plus_so_count.items():
        print(f"{count}: {v} 数量{n}")
        count +=1

if __name__ == "__main__":

    lambda_sp = vector(QQ,[3/2, 3/2])
    lambda_so = vector(QQ,[1/2])
    lambda_sp_plus_so = vector(QQ, [1/2,1/2,1/2])
    lambda_sp_plus_so_next = vector(QQ,[-1/2,1/2,1/2])
    
    n=3 
    m=1
    while True:
        print("*******开始计算***********************************")
        print("你想要做什么？")
        print("1,直接计算")
        print("2,计算K_L分解")
        print("3,知道P,计算P_tensor_W并投射到块lambda")
        print("4,计算集合极小权")
        print("5,计算两个向量集合相减")
        print("6,直接计算typical权")
        print("7,计算P_tensor_W,但是不投射,找到极小权")
        print("8,整理权集")
        print("9,批量直接计算atypical权集的特征标")
        print("10,慎用,一条龙计算,出发点是一个front文件, 一个behind文件, 一个lam文件")

        while True:
            try:
                select_case = int(input("你的选择是: "))
                break  # 如果成功转换为整数，跳出循环
            except ValueError:
                print("输入无效，请输入一个整数。")

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

            which_mod = int(lines[2])
            if which_mod ==1:
                print("使用模:V")# 将输入转换为有理数列表并创建向量
            elif which_mod==2:
                print("使用模:S2V")# 将输入转换为有理数列表并创建向量
            elif which_mod==3:
                print("使用模:g")# 将输入转换为有理数列表并创建向量
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
            P_mu_tensor_V_after_Pr = again_calc(again_lam,P_mu_tensor_V_after_Pr,which_mod,n,m)
        elif select_case==4:
            user_input = input("请输入权集合set所在文档的名字:")# 将输入转换为有理数列表并创建向量
            weight_set = read_vectors_from_file("test//"+user_input+".txt")
            lowest_module = Lowest_Module(n,m)
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
            
            user_input_W = input("请输入模1，2，3:  ")# 将输入转换为有理数列表并创建向量

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
            P_weights = vectors_set_min(front,behind)

#            with open("test://P.txt", "w", encoding="utf-8") as f:
#                f.write("\n".join(P_weights))  # 每个元素一行

#            test = read_vectors_from_file("test//"+"P"+".txt")
#            print(f"写入文件P共{len(test)}个向量")

            with open('test//'+user_input_lam +'.txt', 'r') as f:
                lines = f.readlines()
            rat_list = [QQ(x.strip()) for x in lines[0].split(',')]
            again_lam = vector(QQ,rat_list)          

            which_mod = int(lines[1])
            if which_mod ==1:
                print("使用模:V")# 将输入转换为有理数列表并创建向量
            elif which_mod==2:
                print("使用模:S2V")# 将输入转换为有理数列表并创建向量
            elif which_mod==3:
                print("使用模:g")# 将输入转换为有理数列表并创建向量
            P_mu_tensor_V_after_Pr = again_calc(again_lam, P_weights, which_mod,n,m)
            if P_mu_tensor_V_after_Pr is None:
                print(f"不是极小权，结束循环")
                continue
            print(f"循环要处理的总数: {len(P_mu_tensor_V_after_Pr)}")

 #           with open("test://front2.txt", "w", encoding="utf-8") as f:
 #               f.write("\n".join(P_mu_tensor_V_after_Pr))  # 每个元素一行

            weight_set = P_mu_tensor_V_after_Pr[:]
            lowest_module = Lowest_Module(n,m)
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
                    break

                if weight_set == []:
                    print(f"---- 判断失效，请选择向量重新判断 ----")
                    break

            

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
