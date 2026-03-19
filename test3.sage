# 构造Coxeter群 of type C3 (对应 sp(6))
W = WeylGroup(['A', 3], prefix="w")
w = W.simple_reflections()  # w[1], w[2], w[3]

# 构造x = w3 * w1 * w2
x = w[3] * w[1] * w[2]

# 构造y = w1 * w2 * w3 * w1 * w2
y = w[3] * w[1] * w[2] * w[3] * w[1] * w[2]

# 最长元素 w0
w0 = W.long_element()

# 计算 xw0 和 yw0
xw0 = x * w0
yw0 = y * w0

# 打印长度以验证合理性
print("Length of xw0:", xw0.length())
print("Length of yw0:", yw0.length())

# 检查 Bruhat order: 是否 yw0 <= xw0 ?
if not yw0.bruhat_le(xw0):
    print("yw0 is not <= xw0 in Bruhat order, so KL polynomial is 0.")
    P = 0
else:
    # 构造KL多项式环
    R.<q> = PolynomialRing(QQ)
    KL = KazhdanLusztigPolynomial(W, q)

    # 注意：Sage 的 KL 函数是 KL.P(u, v) = P_{u,v}(q)
    P_poly = KL.P(yw0, xw0)
    print("Kazhdan-Lusztig polynomial P_{yw0, xw0}(q) =", P_poly)
    P = P_poly.subs(q=1)

print("P_{yw0, xw0}(1) =", P)


for w in W:
    if not w.to_matrix() == w.matrix():
        print(f"==={w}")
