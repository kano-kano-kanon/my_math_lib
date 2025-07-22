#my_math_lib_main_ver0.01.py

#定数定義
class const:
    #数学定数
    pi = 3.1415926535897932384                  # 円周率（20桁）
    e = 2.7182818284590452353                   # ネイピア数（20桁）
    phi = 0.6180339887498948482                 # 黄金数（20桁）
    apery = 1.2020569031595942853                # アペリーの定数（20桁）
    gelf = 1.6325269194381528447                 # ゲルフォント・シュナイダー定数（20桁）
    pla = 1.3247179572447460259                  # プラスチック数（19桁 → 真値）
    euler = 0.5772156649015328606                # オイラー定数（20桁）
    catalan = 0.91596559417721901505             # カタラン定数（20桁）
    khinchin = 2.68545200106530644530            # キンチン定数（20桁）
    feigenbaum_d = 4.66920160910299067185        # Feigenbaum定数 δ（20桁）
    feigenbaum_a = 2.50290787509589282228        # Feigenbaum定数 α（20桁）
    glaisher = 1.28242712910062263687             # Glaisher-Kinkelin定数
    euler_gamma = 0.57721566490153286060          # オイラー・マスケローニ定数（重複注意）
    meissel_mertens = 0.26149721284764278375      # Meissel–Mertens定数
    twin_prime_const = 0.66016181584686957392     # 双子素数定数
    brass = 1.45607494858268967145                 # ブラス定数
    yang_lee_edge = 0.29559774252203960404         # Yang–Leeのエッジ定数
    gauss_kuzmin = 0.30366300289873263056          # Gauss-Kuzmin定数
    vilars_constant = 0.76422365358922066299       # Vilarsの定数

    #物理定数（2024 CODATA）
    c = 299792458                                 # 光速 [m/s]（定義値）
    G = 6.67430e-11                               # 万有引力定数 [m^3 kg^-1 s^-2]
    h = 6.62607015e-34                            # プランク定数 [J·s]
    hbar = 1.054571817e-34                         # ディラック定数 ħ = h / (2π)
    k = 1.380649e-23                              # ボルツマン定数 [J/K]
    NA = 6.02214076e23                            # アボガドロ定数 [/mol]
    qe = 1.602176634e-19                          # 素電荷 [C]
    eps0 = 8.8541878128e-12                       # 真空の誘電率 [F/m]
    mu0 = 1.25663706212e-6                        # 真空の透磁率 [N/A²]
    R = 8.314462618                               # 気体定数 [J/mol·K]
    alpha = 0.0072973525693                        # 微細構造定数（無次元）
    Rydberg = 10973731.568160                      # リュードベリ定数 [m^-1]
    mu_B = 9.274009994e-24                         # ボーア磁子 [J/T]
    mu_N = 5.0507837461e-27                        # 核磁子 [J/T]
    eV = 1.602176634e-19                           # 電子ボルト（素電荷と同じ）
    fine_structure = 7.2973525693e-3               # 微細構造定数（無次元、別名）
    hartree_energy = 4.3597447222071e-18           # ハートリーエネルギー [J]
    rydberg_energy = 2.1798723611035e-18            # リュードベリエネルギー [J]
    planck_length = 1.616255e-35                    # プランク長 [m]
    planck_time = 5.391247e-44                       # プランク時間 [s]
    stefan_boltzmann = 5.670374419e-8               # ステファン・ボルツマン定数 [W m^-2 K^-4]
    gas_constant = 8.314462618                       # 気体定数（Rと同値）
    bohr_radius = 5.29177210903e-11                  # ボーア半径 [m]

    #化学定数
    atm = 101325.0                                  # 標準大気圧 [Pa]
    molar_gas_const = 8.314462618                   # 気体定数（Rと同値）
    Faraday = 96485.33212                           # ファラデー定数 [C/mol]
    standard_pressure = 101325                       # 標準圧力 [Pa]
    avogadro = 6.02214076e23                         # アボガドロ定数（NAと同じ）
    boltzmann = 1.380649e-23                         # ボルツマン定数（kと同じ）
# クラス外に全部展開
pi = const.pi
e = const.e
phi = const.phi
apery = const.apery
gelf = const.gelf
pla = const.pla
euler = const.euler
catalan = const.catalan
khinchin = const.khinchin
feigenbaum_d = const.feigenbaum_d
feigenbaum_a = const.feigenbaum_a
glaisher = const.glaisher
euler_gamma = const.euler_gamma
meissel_mertens = const.meissel_mertens
twin_prime_const = const.twin_prime_const
brass = const.brass
yang_lee_edge = const.yang_lee_edge
gauss_kuzmin = const.gauss_kuzmin
vilars_constant = const.vilars_constant
c = const.c
G = const.G
h = const.h
hbar = const.hbar
k = const.k
NA = const.NA
qe = const.qe
eps0 = const.eps0
mu0 = const.mu0
R = const.R
alpha = const.alpha
Rydberg = const.Rydberg
mu_B = const.mu_B
mu_N = const.mu_N
eV = const.eV
fine_structure = const.fine_structure
hartree_energy = const.hartree_energy
rydberg_energy = const.rydberg_energy
planck_length = const.planck_length
planck_time = const.planck_time
stefan_boltzmann = const.stefan_boltzmann
gas_constant = const.gas_constant
bohr_radius = const.bohr_radius
atm = const.atm
molar_gas_const = const.molar_gas_const
Faraday = const.Faraday
standard_pressure = const.standard_pressure
avogadro = const.avogadro
boltzmann = const.boltzmann



def sin(x):
    # xはラジアン
    # 10項までのテイラー展開
    s = 0
    num = x
    sign = 1
    fact = 1
    for i in range(1, 20, 2):
        s += sign * num / fact
        num *= x * x
        fact *= (i + 1) * (i + 2)
        sign *= -1
    return s

def sqrt(x):
    # ニュートン法
    if x < 0:
        raise ValueError("sqrt of negative number")
    if x == 0:
        return 0
    guess = x
    for _ in range(20):
        guess = (guess + x / guess) / 2
    return guess

def cos(x):
    # x in radians, 10 terms Taylor expansion
    s = 0
    num = 1.0
    sign = 1
    fact = 1
    for i in range(0, 20, 2):
        s += sign * num / fact
        num *= x * x
        fact *= (i + 1) * (i + 2)
        sign *= -1
    return s

def tan(x):
    c = cos(x)
    if abs(c) < 1e-12:
        raise ValueError("tan(x): cos(x) is too close to zero")
    return sin(x) / c

def ln(x):
    """高精度自然対数（ビルトインのみ, 20桁目標）"""
    if x <= 0:
        raise ValueError("ln(x) is undefined for x <= 0")
    # x = m * 2^k に変形
    k = 0
    while x > 2:
        x /= 2
        k += 1
    while x < 1:
        x *= 2
        k -= 1
    # ln(x) = k*ln(2) + ln(y) (1 <= y < 2)
    # ln(2)をメリカトリ級数で高精度に求める
    def ln2_mer():
        t = 1  # ln(2) = log(1+1)
        s = 0
        sign = 1
        for n in range(1, 100):
            term = sign * (t**n) / n
            s += term
            sign *= -1
            if abs(term) < 1e-21:
                break
        return s
    ln2 = ln2_mer()
    # x = 1 + t (0 < t < 1)
    t = x - 1
    s = 0
    sign = 1
    for n in range(1, 200):
        term = sign * (t**n) / n
        s += term
        sign *= -1
        if abs(term) < 1e-21:
            break
    return k * ln2 + s

def fact(n):
    res = 1
    for i in range(2, int(n)+1):
        res *= i
    return res

def C(n, r):
    n = int(n)
    r = int(r)
    if r < 0 or n < 0 or r > n:
        return 0
    return fact(n) // (fact(r) * fact(n - r))

def P(n, r):
    n = int(n)
    r = int(r)
    if r < 0 or n < 0 or r > n:
        return 0
    return fact(n) // fact(n - r)

def mean(lst):return sum(lst) / len(lst) if len(lst) > 0 else 0

def deviation(lst):
    m = mean(lst)
    return [x - m for x in lst]

def stddev(lst):
    m = mean(lst)
    if len(lst) < 2:
        return 0
    return (sum((x - m)**2 for x in lst) / (len(lst)-1)) ** 0.5

def sumall(lst):return sum(lst)

def prodall(lst):
    res = 1
    for x in lst:
        res *= x
    return res
