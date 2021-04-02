import copy
import math
import random


#class below
class Matrix():
    def __init__(self, list_2D):
        self.matrix = self.check_and_to_float(list_2D)
        self.row = len(self.matrix)
        self.col = len(self.matrix[0])
        self.shape = (self.row, self.col)
           
    def check_and_to_float(self, list_2D):
        temp = copy.deepcopy(list_2D)
        count = len(temp[0])
        for i in range(0, len(list_2D)):
            if (len(temp[i]) != count):
                raise Exception('输入列表无法构成矩阵，检查每列是否长度相等')
            for j in range(0, len(list_2D[0])):
                temp[i][j] = float(list_2D[i][j])
        return temp
        
    def show(self):
        return copy.deepcopy(self.matrix)#show返回一个拷贝值

    
    def inverse(self):
        '''假设此矩阵为可逆方阵'''
        re = []
        for i in range(0, self.row):
            vec = [0 for j in range(0, self.row)]
            vec[i] = 1
            result = (self / Vector(vec)).show()
            re.append(result)
        return Matrix(matrix_transposition(re))


    def __mul__(self, rhs):
        if (self.col != rhs.row):
            raise Exception('两个矩阵行列不匹配，无法相乘')
        temp1 = self.show()
        temp2 = rhs.show()
        result = [[sum([temp1[i][k]*temp2[k][j] for k in range(0, self.col) ]) 
                        for j in range(0, rhs.col)] for i in range(0, self.row)]
        return Matrix(result)



    #def __len__(self): 
    #len函数无法实现，因为必须返回一个整数值，使用无法同时返回行数和列数
    #所以使用shape属性

    def __add__(self, rhs):
        temp1 = self.show()
        temp2 = rhs.show()
        for i in range(0, self.row):
            for j in range(0, self.col):
                temp1[i][j] += temp2[i][j]
        return Matrix(temp1)

    def __sub__(self, rhs):
        temp1 = self.show()
        temp2 = rhs.show()
        for i in range(0, self.row):
            for j in range(0, self.col):
                temp1[i][j] -= temp2[i][j]
        return Matrix(temp1)

    def __truediv__(self, rhs):
        if (self.row != self.col):
            raise Exception('系数矩阵必须是方阵')
        if (self.row != len(rhs)):
            raise Exception('系数矩阵和向量必须具有相同的规模')
        results =  solve_equations_by_gauss(self.show(), rhs.show())
        if results == -1:
            return '方程组没有唯一解'
        else:
            return Vector(results)

    def __repr__(self):
        return 'Matrix(' + str(self.show()) + ')'

    def simple_solve_equations_by_gauss(self, rhs):
        if (self.row != self.col):
            raise Exception('系数矩阵必须是方阵')
        if (self.row != len(rhs)):
            raise Exception('系数矩阵和向量必须具有相同的规模')
        results =  simple_solve_equations_by_gauss(self.show(), rhs.show())
        if results == -1:
            return '对角线上有零，此法不通，必须选主元'
        else:
            return Vector(results)

    def solve_equations_by_cholesky(self, rhs):
        if (self.row != self.col):
            raise Exception('系数矩阵必须是方阵')
        if (self.row != len(rhs)):
            raise Exception('系数矩阵和向量必须具有相同的规模')
        results =  solve_equations_by_cholesky(self.show(), rhs.show())
        return Vector(results)

    def improved_solve_equations_by_cholesky(self, rhs):
        if (self.row != self.col):
            raise Exception('系数矩阵必须是方阵')
        if (self.row != len(rhs)):
            raise Exception('系数矩阵和向量必须具有相同的规模')
        results =  improved_solve_equations_by_cholesky(self.show(), rhs.show())
        return Vector(results)

    def solve_LS_problem_by_QR(self, rhs):
        if (self.row != len(rhs)):
            raise Exception('系数矩阵行数和向量元素个数必须相同')
        results = solve_LS_problem_by_QR(self.show(), rhs.show())
        return Vector(results)


class Vector():
    def __init__(self, list_1D):
        self.vector = [float(i) for i in list_1D]
        self.size = len(list_1D)
    
    def show(self):
        return self.vector[:]

    def inf_norm(self):
        return max([abs(i) for i in self.vector])

    def __len__(self):
        return self.size
    
    def __add__(self, rhs):
        temp1 = self.show()
        temp2 = rhs.show()
        for i in range(0, self.size):
                temp1[i] += temp2[i]
        return Vector(temp1)

    def __sub__(self, rhs):
        temp1 = self.show()
        temp2 = rhs.show()
        for i in range(0, self.size):
                temp1[i] -= temp2[i]
        return Vector(temp1)
    
    def __repr__(self):
        return 'Vector(' + str(self.show()) + ')'

'''
functions below are provided to classes,
in functions below, i assume that the inputs are behaved and
any exception is preprocessed in the methods of class.
'''

def solve_equations_by_gauss(A, b):
    '''A为n阶方阵，
       b为n维向量，
       A为奇异时返回-1'''
    results = gauss_elimination(A, b)
    if (results != -1):
        A, b = results
    else:
        return -1
    #'unit'表示这是单位下三角阵，目的在于不必把A拆成L和U
    y = solve_lower_equations(A, b, 'unit')
    x = solve_upper_equations(A, y)
    return x

def solve_equations_by_cholesky(A, b):
    '''A为n阶对称正定阵，
       返回解x'''
    L = get_cholesky_factor(A)
    y = solve_lower_equations(L, b)
    x = solve_upper_equations(matrix_transposition(L), y)
    return x

def improved_solve_equations_by_cholesky(A, b):
    '''A为n阶对称正定阵，
       返回解x'''
    n = len(A)
    L = improved_cholesky(A)
    z = solve_lower_equations(L, b, 'unit')
    y = [z[i] / A[i][i] for i in range(0, n)]
    #不再需要D了，将L对角线上元素置为1，方便求上三角解
    for i in range(0, n):
        A[i][i] = 1
    x = solve_upper_equations(matrix_transposition(L), y)
    return x

def solve_LS_problem_by_QR(A, b):
    ''''''
    re = QR_elimination(A)
    R = re[3]
    Q = re[2]
    Q1_T = matrix_transposition([Q[i][0:len(R)] for i in range(0, len(Q[0]))])
    c1 = [sum([Q1_T[i][j] * b[j] for j in range(0, len(b))]) for i in range(0, len(Q1_T))]
    x = solve_upper_equations(R, c1)
    #print(R,'\n', Q,'\n', Q1_T,'\n', c1,'\n', x)
    #print(x)
    return x

def matrix_transposition(A):
    m, n = len(A), len(A[0])
    B = []
    for i in range(0 ,n):
        col = [A[k][i] for k in range(0, m)]
        B.append(col)
    return B

def solve_lower_equations(A, b, sign='general'):
    '''A为n阶下三角非奇异方阵，
       b为n维向量，
       sign='unit'表示这是单位下三角阵, 默认sign='general'表示这是一般下三角阵
       此法用于解决Lx=b'''
    n = len(A)
    if (sign == 'general'):
        for j in range(0, n):
            if (A[j][j] != 0):
                b[j] = b[j] / A[j][j]
            else:
                return -1
            for i in range(j+1, n):#取到n-1，当j+1=n时，此循环消失
                b[i] = b[i] - b[j] * A[i][j]
        return b
    else:
        #当A为单位下三角阵时
        for j in range(0, n):
            for i in range(j+1, n):#取到n-1，当j+1=n时，此循环消失
                b[i] = b[i] - b[j] * A[i][j]
        return b


def solve_upper_equations(A, b):
    '''B为n阶上三角非奇异方阵,
       b为n维向量,
　　　 此法用于解决Ux=b'''
    n = len(A)
    for j in range(n-1, -1, -1):
        if (A[j][j] != 0):
            b[j] = b[j] / A[j][j]
        else:
            return -1
        for i in range(j-1, -1, -1):
            b[i] = b[i] - b[j] * A[i][j]
    return b

def gauss_elimination(A, b=0):
    '''A为n阶方阵，
       b为n维向量或者不输入，
       能力所限，A为奇异时返回-1'''
    n = len(A)
    for k in range(0, n-1):#第n-1列不必再管
        max_element = k
        for p in range(k+1, n):
            if (abs(A[p][k]) > abs(A[max_element][k])):
                max_element = p
        #已确定最大行开始交换
        A[k], A[max_element] = A[max_element], A[k] #交换列表引用
        if (b != 0):
            #在gauss分解时不记录置换矩阵，直接对b进行置换
            b[k], b[max_element] = b[max_element], b[k]
        if (A[k][k] != 0.0):
            for i in range(k+1, n):
                A[i][k] = A[i][k] / A[k][k]
                for j in range(k+1, n):
                    A[i][j] = A[i][j] - A[i][k] * A[k][j]
        else:
            return -1
    if (A[n-1][n-1] == 0.0):
        return -1
    if (b != 0):
        return [A, b]
    else:
        return A

def simple_gauss_elimination(A):
    '''A为n阶方阵，
       遇到除0返回-1'''
    n = len(A)
    for k in range(0, n-1):#第n-1列不必再管
        if (A[k][k] != 0.0):
            for i in range(k+1, n):
                A[i][k] = A[i][k] / A[k][k]
                for j in range(k+1, n):
                    A[i][j] = A[i][j] - A[i][k] * A[k][j]
        else:
            return -1
    if (A[n-1][n-1] == 0.0):
        return -1
    else:
        return A

def simple_solve_equations_by_gauss(A, b):
    '''A为n阶方阵，
       b为n维向量，
       简单LU分解除0返回-1'''
    results = simple_gauss_elimination(A)
    if (results != -1):
        A = results
    else:
        return -1
    #'unit'表示这是单位下三角阵，目的在于不必把A拆成L和U
    y = solve_lower_equations(A, b, 'unit')
    x = solve_upper_equations(A, y)
    return x

def get_cholesky_factor(A):
    '''A为n阶对称正定阵，
       返回矩阵L'''
    n = len(A)
    for k in range(0, n):
        A[k][k] = math.sqrt(A[k][k])
        for i in range(k+1, n):
            A[i][k] = A[i][k] / A[k][k]
        for j in range(k+1, n):
            for i in range(j, n):
                A[i][j] = A[i][j] - A[j][k]*A[i][k]
    return A
        
def improved_cholesky(A):
    '''A为n阶对称正定阵，
       返回L和D的叠加，合在一个阵中'''
    n = len(A)
    for j in range(0, n):
        v = [A[i][i]*A[j][i] for i in range(0, j)]
        for i in range(0, j):
            A[j][j] -= A[j][i]*v[i]
        for i in range(j+1, n):
            for k in range(0, j):
                A[i][j] -= A[i][k]*v[k]
            A[i][j] /= A[j][j]
    return A


def unit_matrix(n):
    re = [[0 for i in range(0, n)] for j in range(0, n)]
    for i in range(0, n):
        re[i][i] = 1
    return re

def house_holder(x):
    '''x为一个一维数组，返回在最后附上beta的v'''
    n = len(x)
    inf_norm = Vector(x).inf_norm()
    x = [x[i] / inf_norm for i in range(0, n)]
    sigma = sum([x[i]*x[i] for i in range(1, n)])
    v = [0 for i in range(0, n)]
    v[1:-1] = [x[i] for i in range(1, n)]
    beta = 0
    if (sigma == 0):
        beta = 0
    else:
        arf = math.sqrt(x[0]*x[0] + sigma)
        if (x[0] <= 0):
            v[0] = x[0] - arf
        else:
            v[0] = -1*sigma/(x[0]+arf)
        beta = 2*v[0]*v[0]/(sigma+v[0]*v[0])
        a = v[0]
        v = [v[i]/a for i in range(0, n)]
    v.append(beta)
    return v

def QR_elimination(A):
    '''一般假设输入的二维数组A行数大于等于列数
       返回A的QR分解，A中存储所以v和R，另外存储与v对应的Beta在d中
       返回[A, d]'''
    m = len(A)
    #不考虑A为空
    n = len(A[0])
    d = []
    H = []#保存待处理的Hk，此处Hk形状还未补全
    for j in range(0, n):#注意这里是n不是n-1，害得我少算了一个搞了半天
        if (j < m-1):#注意这里是m-1而不是m
            temp1 = house_holder([A[i][j] for i in range(j, m)])
            v, beta = temp1[0:-1], temp1[-1]
            process_matrix = Matrix([A[i][j:n] for i in range(j, m)])
            temp2 = Matrix(unit_matrix(m-j))
            temp3 = Matrix([[beta*v[k]*v[i] for i in range(0, m-j)] for k in range(0, m-j)])
            H.append((temp2-temp3).show());
            result = ((temp2 - temp3)*process_matrix).show()
            for i in range(j, m):
                A[i][j:n] = result[i-j]
            d.append(beta)
            for i in range(j+1, m):
                A[i][j] = v[i-j]
    #为了方便我们这里花大量时间把具体的Q算出来
    R = copy.deepcopy(A)
    R = [R[i] for i in range(0, len(R[0]))]
    for i in range(0, len(R)):
        for j in range(0, i):
            R[i][j] = 0
    process_num = len(d)
    if (process_num == 0):
        return [A, d, [], R]
    #下面处理H，把形状补全成m*m，得到H1,H2...
    H = [Matrix(process_H(H[i], m)) for i in range(0, len(H))]
    if (len(H) == 1):
        Q = H[0].show()
    else:
        H_left = H[0]
        for i in range(1, len(H)):
            H_left = H_left * H[i]
        Q = H_left.show()
    return [A, d, Q, R]


def iterative_method(M, g, x0, max_num):
    '''M为n阶迭代矩阵(方阵)，g为常数项，x0为初始向量，max_num为最大迭代次数'''
    n = len(M)
    x1 = [k for k in x0]
    x2 = []
    for i in range(0, max_num):
        x2 = [sum([M[i][j] * x1[j] for j in range(0, n)])+g[i] for i in range(0, n)]
        temp = max([abs(x1[i] - x2[i]) for i in range(0, n)])

        if (temp <= 0.00001 ):
            return [x2, i+1]
        else:
            x1 = x2
    return [x2, max_num]    

def conjugate_gradient_method(A, x, b, km):
    '''A为对称正定阵，x为初始向量，b为常向量，km为最大迭代次数，返回方程组Ax=b的解x和迭代次数k'''
    k = 0
    n = len(b)
    r = [b[j] - sum([A[j][i]*x[i] for i in range(0, n)]) for j in range(0, n)]
    q_new = sum([r[i]*r[i] for i in range(0, n)])
    q_old = 0
    p = []
    while ((math.sqrt(q_new) > 0.000001) and k < km):
        k += 1;
        if (k == 1):
            p = r[:]
        else:
            beta = q_new / q_old
            p = [r[i] + beta*p[i] for i in range(0, n)]
        w = [sum([A[j][i]*p[i] for i in range(0, n)]) for j in range(0, n)]
        arf = q_new/sum([p[i]*w[i] for i in range(0, n)])
        x = [x[i]+arf*p[i] for i in range(0, n)]
        r = [r[i]-arf*w[i] for i in range(0, n)]
        q_old = q_new
        q_new = sum([r[i]*r[i] for i in range(0, n)])
    return [x, k]

def max_module_root_of_polynomial(coe):
    '''输入一个系数向量coe，coe[0]表示x^(n-1)的系数，一直到常数项，
    x^n系数默认为1，输出模最大根'''
    n = len(coe)
    A = [[0 for i in range(0, n)] for j in range(0, n)]
    A[0][0] = -1*coe[0]
    for i in range(1, n):
        A[i][i-1] = 1
        A[0][i] = -1*coe[i]
    re = power_method(A)
    return re[0]

def power_method(A, num=1000):
    '''输入方阵A，返回矩阵A的模最大特征值和特征向量'''
    n = len(A)
    u = [random.random() for i in range(0, n)]
    u[0] = 1
    mu = 0
    k = 1
    while (k <= num):
        y = [sum([A[i][j] * u[j] for j in range(0, n)]) for i in range(0, n)]
        mu = max([abs(y[i]) for i in range(0, n)])
        u = [y[i]/mu for i in range(0, n)]
        k += 1
    return [mu, u]


"""
Functions below are independent, providing extra effect.
"""

def special_norm_1_by_optimization(A):
    '''A为n阶矩阵，本方法返回A^(-T)的1范数，而不是A的1范数'''    
    k = 1
    n = len(A)
    x = [1/n for i in range(0, n)]
    while (k == 1):
        w = (Matrix(matrix_transposition(A)) / Vector(x)).show()
        v = [1 if i>0 else 0 if i==0 else -1  for i in w]
        z = (Matrix(A) / Vector(v)).show()
        if max([abs(i) for i in z]) <= sum([z[i]*x[i] for i in range(0, len(z))]):
            v = sum([abs(i) for i in w])
            k = 0
        else:
            x = [0 for i in range(0, n)]
            j = 0
            for i in range(0, n):
                if (abs(z[j]) < abs(z[i])):
                    j = i
            x[j] = 1
            k = 1
    return v

def process_H(H, m):
    '''这个函数仅给QR_elimination(A)用'''
    re = [[0 for i in range(0, m)] for j in range(0, m)]
    for i in range(0, m-len(H)):
        re[i][i] = 1
    for i in range(m-len(H), m):
        re[i][m-len(H):] = H[i-m+len(H)]
    return re
        

if __name__ == '__main__':
    print('选主元高斯消元法Example:')
    A = Matrix([[3, 5, 7], [2, 4, 8], [1, 7, 5]])
    b = Vector([1, 3, 9])
    print('A = ', A)
    print('b = ', b)
    print('Ax = b 解得 x = ', A / b)
    print('\n')
    
    print('对称正定阵的平方根法Example:')
    A = Matrix([[16, 4, 8, 4], [4, 10, 8, 4], [8, 8, 12, 10], [4, 4, 10, 12]])
    b = Vector([32, 26, 38, 30])
    print('A = ', A)
    print('b = ', b)
    print('Ax = b 解得 x = ', A.improved_solve_equations_by_cholesky(b))
