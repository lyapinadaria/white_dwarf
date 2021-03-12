# -*- coding: utf-8 -*-

from __future__ import unicode_literals

import numpy as np
import matplotlib.pyplot as plt
from math import exp
import random

#  ФУНКЦИЯ ПОДСЧЁТА СТРОК В ФАЙЛЕ
def count_lines(filename, chunk_size=1<<13):
    """ Функция считает строки в файле 
    
    filename : имя файла
    
    Возвращает количество строк в файле
    """
    with open(filename) as file:
        return sum(chunk.count('\n')
                   for chunk in iter(lambda: file.read(chunk_size), ''))
   
    

#  ФУНКЦИЯ ВЫДАЧИ ИМЕНИ ФАЙЛА С ДАННЫМИ
def get_m_file(b):
    """ Функция выдаёт имя файла по полученному значению массы
    
    b : параметр массы
    
    Возвращает имя файла с данными
    """
    return {
                b < 0.04:  'ResM04.txt',
        0.04 <= b < 0.15:  'ResM05.txt',
        0.15 <= b < 0.52:  'ResM06.txt',
        0.52 <= b < 0.81:  'ResM07.txt',
        0.81 <= b < 0.91:  'ResM08.txt',
        0.91 <= b < 0.94:  'ResM09.txt',
        0.94 <= b:         'ResM10.txt' 
       }[True]


#  ФУНКЦИЯ ОПРЕДЕЛЕНИЯ ПОГЛОЩЕНИЯ
def get_Nh(R0):
    """ Функция определяет поглощение в зависимости от расстояния в модели Local Bubble
    
    R : расстояние 
    
    Возвращает число атомов на луче зрения
    """
    PC=3.086e18 # Парсек
    LB_Line=100 # Граница Local Bubble
    Nh_coef_in=0.1
    Nh_coef_out=1.
    
    if R0 < LB_Line * PC:
        Nh = Nh_coef_in * R0 
    else:
        Nh=Nh_coef_in * LB_Line * PC  + Nh_coef_out * (R0 - LB_Line * PC) 
    return Nh

    
# ФУНКЦИЯ ОПРЕДЕЛЕНИЯ ИНДЕКСА (Служебная для аппроксимации)
def get_index(Age,e):
    """ Функция определяет индекс в списке возрастов
    
    Age : список возрастов
    e : возраст, после которого шаг меняется
    
    Возвращает номер последней строки перед сменой шага
    """
    for i in enumerate(Age): 
        if i[1] > e:  
            return i[0]-1    

# ФУНКЦИЯ ОПРОКСИМАЦИИ (Служебная для аппроксимации)

def aprox(Y,X,index,X3):
    """ Функция аппроксимации. Проводит линейную аппроксимацию.
    
    Y: массив температур
    X : массив возрастов
    index : ближайшая строка в массиве X и Y через которую мы считаем аппроксимацию 
    X3 : возраст, для которого нужно посчитать температуру 
    
    Возвращает температуру для возраста в Х3
    """
    T_Teff = (Y[index+1] - Y[index]) /  (X[index+1] - X[index]) * (X3 - X[index]) +  Y[index]
    return T_Teff   

# ФУНКЦИЯ ВЫДАЧИ ИНДЕКСА ДЛЯ КОЭФФИЦИЕНТОВ
def get_c_i(E):
    """ Функция выдачи индекса для коэффициентов для поглощения в зависимости от энергии C0C1C2
    
    Е: энергия
  
    
    Возвращает номер строки в списке коэффициентов
    """
    return {
                 E < 0.284:  0,
        0.284 <= E < 0.400:  1,
        0.400 <= E < 0.532:  2,
        0.532 <= E < 0.707:  3,
        0.707 <= E < 0.867:  4,
        0.867 <= E < 1.303:  5,
        1.303 <= E < 1.840:  6,
        1.840 <= E < 2.471:  7,
        2.471 <= E < 3.210:  8,
        3.210 <= E < 4.038:  9,
        4.038 <= E < 7.111: 10,
        7.111 <= E < 8.331: 11, 
        8.331 <= E:         12
       }[True]


# ФУНКЦИЯ ПОДСЧЁТА LogN / LogS
def get_logNS(N,Nph,q,Smax,EXPOS):
    """ Функция подсчёта LogN / LogS
    
    N: Число фотонов от объекта
    Nph: Поток фотонов от объекта
    q: Массив числа объектов
    Smax: Поток фотонов
    EXPOS: Экспозиция
  
    Функция заполняет массив q
    """
    for x in range(Smax): 
        if N * EXPOS > Nph[x]:
            q[x] = q[x] + 1
        else:
            N=0
            break      
   
# ФУНКЦИЯ ПОСТРОЕНИЯ ГРАФИКА
def get_graf():
    """ Функция построения графика
    
  
    Функция строит график
    """
    Nph,q = np.genfromtxt('N(S).txt',unpack=True)
    plt.figure(1)
    plt.figure(figsize=(12, 12), dpi=400)
    plt.xlabel('N')
    plt.ylabel('Q')
    plt.grid(True)  
    plt.yscale('log')
    plt.xscale('log')
    plt.plot(Nph, q, 'v-.g', mec='r')
    plt.show()



def glamurny_carlik():
    """ Функция построения популяционного синтеза. Функция синтезирует популяцию... зачем... 
    
    Функция создаёт файл 'N(S).txt'. Функция строит график N(S).
    где: N- число объектов
         S- поток фотонов
    """
    RC=count_lines('radius.txt')
    ES=count_lines('S+sigma.txt')
    print('Start WhiteDwarf')
    print('RC=',RC)
    print('ES=',ES)
    #КОНСТАНТЫ
    Smax=10000 #Smax=620000 Поток фотонов
    WDCount=int(1e3) #WDCount=10000 Количество белых карликов 
    h=1.05e-27
    c=3.e10
    k=1.38e-16   
    Em,S=np.genfromtxt('S+sigma.txt',unpack=True) #выгружаем эффективную площадь
    kev_to_erg = 1.6e-9
    R=np.genfromtxt('radius.txt',unpack=True) # распределение по расстояниям
    N=0 #число фотонов для данного объекта
    Nph = np.arange(10, Smax * 10 + 1, 10) #массив значений энергии для logNlogS
    q=np.zeros(Smax)#массив для числа объектов

    EXPOS = 2500 #Экспозиция

    C0C1C2 = [[34.6,    267.9, -476.1], #0
          [78.1,    18.8,     4.3], #1
          [71.4,    66.8,   -51.4], #2
          [95.5,    145.8,  -61.1], #3
          [308.91, -380.6,   294.], #4
          [120.6,   169.3,  -47.7], #5
          [141.3,   146.8,  -31.5], #6
          [202.7,   104.7,    -17], #7
          [342.7,    18.7,      0], #8
          [352.2,    18.7,      0], #9
          [433.9,    -2.4,   0.75], #10
          [629,      30.9,      0], #11
          [701.2,    25.2,      0]] #12

    
    for i in range(WDCount):
        #print(' ')
        #print('Итерация',i)
        R0 = R[random.randint(0,RC-1)] # Выбираем случайную строку из файла R (расстояние до объекта)
        
        #  ОПРЕДЕЛЯЕМ ПОГЛОЩЕНИЕ
        Nh=get_Nh(R0)
        #print('Nh =',Nh)  # Вывод 
        #print('R0 =',R0)  # Вывод расстояния до объекта
        
        #  ОПРЕДЕЛЯЕМ МАССУ
        b=random.random() # Вероятность попадения в определененный диапазон
        #e=random.randint(10000,100000000) # определяем диапазон возрастов
        Age,r,Teff,L = np.genfromtxt(get_m_file(b),unpack=True)
        #print('m =',b) 
        #print('File =',get_m_file(b))    
        
        #  ОПРЕДЕЛЯЕМ ТЕМПЕРАТУРУ, РАДИУС, СВЕТИМОСТЬ
        e=random.randint(Age[0],Age[-1])  # Выбираем возраст (Границы возраста берем из файла)
        a = get_index(Age,e) # Находим индекс ближайшего элемента слева
        T_Teff = aprox(Teff,Age,a,e) # Аппроксимируем температуру объекта
        r_r = aprox(r,Age,a,e)*1e4 # Аппроксимируем радиус объекта и переводим в [см]
        L_L = aprox(L,Age,a,e) # Аппроксимируем светимость объекта
        #print('teff =',T_Teff)   
        #print('r =',r_r)    
        #print('L =',L_L)
    
    
        # ИНТЕГРАЛ
        for j in range(ES) :  #вылетает из-за слишком большой экспоненты  
             
            # ОПРЕДЕЛЯЕМ КОЭФФИЦИЕНТЫ ДЛЯ ПОГЛОЩЕНИЯ
            c0 = C0C1C2[get_c_i(Em[j])][0] 
            c1 = C0C1C2[get_c_i(Em[j])][1]
            c2 = C0C1C2[get_c_i(Em[j])][2]
          
            # ОПРЕДЕЛЯЕМ ПОГЛОЩЕНИЕ
            sigma = (c2 / Em[j] + c1 / (Em[j] * Em[j]) + c0 / (Em[j] * Em[j] * Em[j])) / 1e24 * Nh 
        
            # ПОДСЧЁТ СУММЫ ДЛЯ ИНТЕГРАЛА 
            if ((Em[j] * kev_to_erg) / (k * T_Teff)) > 680: # Выход из цикла "ИНТЕГРАЛ" 
                 break
            else:
                    N = N + ( 2 * r_r**2 * (Em[j] * kev_to_erg)**2 / R0**2 / h**3 / c**2 / (exp((Em[j] * kev_to_erg) / k / T_Teff) - 1) ) * exp(-sigma) * S[j] * (Em[j+1] - Em[j]) * kev_to_erg #считаем число фотонов
    
    
        #print('N=',N * EXPOS)
        
        # ПОДСЧЁТ LogN / LogS
        get_logNS(N,Nph,q,Smax,EXPOS)
        
        N=0 # Защита на случай переполнения диапазона по потокам
    
    
    # ВЫВОД ИНФОРМАЦИИ В ФАЙЛ    
    
    with open('N(S).txt','w',encoding='utf-8') as f:
       for z in range(Smax-1): #вывод в файл
        f.write(str(Nph[z])+' ')
        f.write(str(q[z])+'\n')
    f.close()

    # СТРОИМ ГРАФИК
    get_graf()
    print('End WhiteDwarf')

