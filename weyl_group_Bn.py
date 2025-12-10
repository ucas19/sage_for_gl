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
class Weyl_Group_Bn:
    """
    Weyl群 W(B_n) = W(so(2n+1)) 的实现，同构于 (Z₂^n) ⋊ S_n。
    群元素可以表示为 (ε, σ)，其中 ε ∈ {±1}^n，σ ∈ S_n。
    群乘法为：(ε1, σ1) * (ε2, σ2) = (ε1 * (σ1(ε2)), σ1 ∘ σ2)
    其中 σ(ε) 表示对 ε 的分量进行置换：σ(ε)_i = ε_{σ^{-1}(i)}
    """   
    def __init__(self,n):
        """
        初始化：n为李代数的秩为n
        """
        self.n = n
        self.W = WeylGroup(['A', n-1], prefix="w")  # Weyl群，生成元命名为 s1, s2, ...
 #       self.Z2_n = AbelianGroup([2]*n, names = "e")
        self.S_n = SymmetricGroup(n)
 #       self.G = SemidirectProductGroup(self.Z2_n, self.S_n) 
        self.W_cox = CoxeterGroup(self.W.cartan_type(), base_ring = AA) 

#    def bn_to_semidirect(self, w):
#    # 获取 w 的矩阵表示
#   #mat = w.to_matrix()
#        v = vector(QQ, range(1,self.n+1))
#        # 提取符号改变信息
#        result = w.to_matrix() * v
##        print(f" {w_inv.to_matrix()} * {v} = {result_inv}")
#        epsilon = []
#        for i in range(self.n):
#            if result[i] > 0:
#                epsilon.append(0)  # 无符号改变
#            else:  # mat[i, i] == -1
#                epsilon.append(1)  # 有符号改变
#        #print(result,epsilon)
#        # 提取置换信息
#        abs_values = [abs(x) for x in result]
#        # 创建置换
#        sigma = self.S_n(abs_values)
##        print(f"{sigma}")
#        return (self.Z2_n(epsilon), sigma.inverse())
#
#    def semidirect_to_bn(self, sd_element):
#        epsilon, sigma = sd_element
#        # 将 AbelianGroup 元素转换为列表
#        eps_list = epsilon.list()
#    #    print(f"测试{eps_list}")
#        #print(sigma)
#        # 创建对角矩阵（符号改变部分）
#        #diag_mat = diagonal_matrix([1 if e == 0 else -1 for e in eps_list])
#        #print((f"eps_list={eps_list}"))
#        # 创建置换矩阵
#        perm_mat = matrix(self.n, self.n, 0)
#        for i in range(self.n):
#            j = sigma(i+1) - 1  # S_n 中的元素作用于 1..n，矩阵索引是 0..n-1
#            perm_mat[j, i] = (-1)**eps_list[j]
#    #        print(f"{j},{i},{perm_mat[j,i]}")
#    #        print(f"{i+1} {j+1}")
#    #    print(f"矩阵是\n{ perm_mat}")
#        # 在 Weyl 群中找到对应的元素
#        for w in self.W:
#            if w.to_matrix() == perm_mat:
#                return w
#        return None  # 如果找不到（理论上不应该发生）
