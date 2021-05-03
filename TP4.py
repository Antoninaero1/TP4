# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 15:03:47 2021

@author: antoc
"""
from math import log #sur python log=ln
from math import exp
from math import sin
from math import cos
import matplotlib.pyplot as plt
import numpy as np


def y0 (x):
    return x**4+3*x-9
def dy0 (x):
    return 4*x**3+3

def y1 (x):
    return x*exp(x)-7
def dy1(x):
    return exp(x)+x*exp(x)

def y2 (x):
    return log(x**2+4)*exp(x)-10
def dy2(x):
    return exp(x)*(log(x**2+4)+((2*x)/(x**2+4)))

def y3(x):
    return 1+sin(x)-2*x
def dy3(x):
    return cos(x)-2

def g0 (x):
    return (9-3*x)**(1/4)
def g1 (x):
    return(log(7/x))
def g2 (x):
    return log(10/log(x**2 +4))
def g3 (x):
    return (1+sin(x))/2





def secant_method(f, x0, x1, epsilon, Nitermax):
    
    n=0
    en=(abs(x1-x0))
    l_n=[n]
    l_xn=[x0]
    l_en=[en]
    
    while en>epsilon and n<=Nitermax :
        n=n+1
        x2=(x0*f(x1)-x1*f(x0))/(f(x1)-f(x0))
        x0=x1
        x1=x2
        en=(abs(x1-x0))
        l_n.append(n)
        l_xn.append(x1)
        l_en.append(en)

    return l_n,l_en,l_xn



def dicho (f,a0,b0,epsilon,Nitermax) :
    
    n=0
    a=a0
    b=b0
    en=(abs(b-a))
    l_n=[n]
    l_xn=[a]
    l_en=[en]

    while en>epsilon and n<=Nitermax :
       m=(a+b)/2
       n+=1
       if f(a)*f(m)<=0 :
           b=m
       else :
           a=m
       en=(abs(b-a))
       l_n.append(n)
       l_xn.append(a)
       l_en.append(en)
       
    return l_n,l_en,l_xn



def ptf(g,X0,epsilon,Nitermax):
    """ Programme du point fixe qui compare l'ecart entre g(xn)[x1] et g(xn+1) [x2] 
    cet écart est comparé à epsilon qui nous permet de définir la précision 
    de la solution g(alpha)=alpha g est deduit d'una fonction f(alpha)=0 
    on inplémente aussi un compteur d'itération pour sorir de la boucle si il n'y a pas de solution trouvé au bout de Nitermax essai
    cela nous permet aussi de savoir à quelle vitesse la méthode du point fixe trouve alpha 
    la fonction return trois listes contenant n l erreur ainsi que l alpha trouvé a chaque itération
    """
    
    n=0
    x1=X0
    en=abs(g(x1)-x1)
    l_n=[n]
    l_xn=[x1]
    l_en=[en]
    
    while en>epsilon and n<Nitermax:
        x2=g(x1)
        n+=1
        en=abs(x2-x1)
        x1=x2
        l_n.append(n)
        l_xn.append(x1)
        l_en.append(en)
        
    return l_n,l_en,l_xn

    

def Newton (f,fder,X0,epsilon,Nitermax):
    """ Programme de la méthode de Newton qui compare l'ecart entre [x1] et g(xn) [x2] 
    Avec g(xn)=x1-(f(x1))/(f'(x1))
    cet écart est comparé à epsilon qui nous permet de définir la précision 
    de la solution alpha tel que f(alpha)=0 
    on inplémente aussi un compteur d'itération pour sorir de la boucle si il n'y a pas de solution trouvé au bout de Nitermax essai
    cela nous permet aussi de savoir à quelle vitesse la méthode de newton trouve alpha
    la fonction return trois listes contenant n l erreur ainsi que l alpha trouvé a chaque itération 
    """
    n=0
    x1=X0
    en=abs(-(f(x1)/fder(x1)))
    l_n=[n]
    l_xn=[x1]
    l_en=[en]    

    while en>epsilon and n<Nitermax:
        x2=x1-(f(x1)/fder(x1))
        n+=1 
        x1=x2
        en=abs(-(f(x1)/fder(x1)))
        l_n.append(n)
        l_xn.append(x1)
        l_en.append(en)
    return l_n,l_xn,l_en
    
data = dicho(y0,0,2,1e-10,5e4)
data1 = secant_method(y0,0,2,1e-10,5e4)
data2 = Newton(y0,dy0,1,1e-10,5e4)
data3 = ptf(g0,0,1e-10,5e4)
plt.semilogy(data[0],data[1])
plt.semilogy(data1[0],data1[1])
plt.semilogy(data2[0],data2[2])
plt.semilogy(data3[0],data3[1])
plt.title('f(x)=x**4+3*x-9')
plt.xlabel('Nombre itération')
plt.ylabel('Erreur')
plt.legend(['Dichotomie','Secante','Newton',"Point fixe"])
plt.show()

data01 = dicho(y1,0,2,1e-10,5e4)
data11 = secant_method(y1,0,2,1e-10,5e4)
data21 = Newton(y1,dy1,1,1e-10,5e4)
data31 = ptf(g1,1.5,1e-10,5e4)
plt.semilogy(data01[0],data01[1])
plt.semilogy(data11[0],data11[1])
plt.semilogy(data21[0],data21[2])
plt.semilogy(data31[0],data31[1])
plt.title('f(x)=x*exp(x)-7')
plt.xlabel('Nombre itération')
plt.ylabel('Erreur')
plt.legend(['Dichotomie','Secante','Newton',"Point fixe"])
plt.show()

data011 = dicho(y2,0,2,1e-10,5e4)
data111 = secant_method(y2,0,2,1e-10,5e4)
data211 = Newton(y2,dy2,1,1e-10,5e4)
data311 = ptf(g2,2,1e-10,5e4)
plt.semilogy(data011[0],data011[1])
plt.semilogy(data111[0],data111[1])
plt.semilogy(data211[0],data211[2])
plt.semilogy(data311[0],data311[1])
plt.title('f(x)=log(x**2+4)*exp(x)-10')
plt.xlabel('Nombre itération')
plt.ylabel('Erreur')
plt.legend(['Dichotomie','Secante','Newton',"Point fixe"])
plt.show()

data0111 = dicho(y3,0,2,1e-10,5e4)
data1111 = secant_method(y3,0,2,1e-10,5e4)
data2111 = Newton(y3,dy3,1,1e-10,5e4)
data3111 = ptf(g3,1,1e-10,5e4)
plt.semilogy(data0111[0],data0111[1])
plt.semilogy(data1111[0],data1111[1])
plt.semilogy(data2111[0],data2111[2])
plt.semilogy(data3111[0],data3111[1])
plt.title('f(x)=1+sin(x)-2*x')
plt.xlabel('Nombre itération')
plt.ylabel('Erreur')
plt.legend(['Dichotomie','Secante','Newton',"Point fixe"])
plt.show()


#print(dicho(y0,0,2,1e-10,5e4),'d1') 

print(Newton(y1,dy0,1,1e-10,5e4),'d2')  
   
#print(dicho(y2,0,2,1e-10,5e4),'d3')    
   
#print(secant_method(y0,0,2,1e-10,5e4),'s1')

#print(secant_method(y1,0,2,1e-10,5e4),'s2')

#print(secant_method(y2,0,2,1e-10,5e4),'s3')
