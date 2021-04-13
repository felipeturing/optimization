using LinearAlgebra #para el producto interno (dot)

function eli_vector(x,D,d)#devuelve el centro x(t+1) dado x(t)=x
    n=length(d)
    return x-(1/(n+1))D*d/(sqrt(transpose(d)*D*d))
end

function eli_matr(x,D,d)#devuelve la matriz D(t+1) dado D(t)=D
    n=length(d)
    return (n^2/(n^2-1))*(D-(2/(n+1))*((D*d)*transpose(D*d))/(transpose(d)*D*d))
end

function pertenece(A,b,x)
    n=length(A[1,:]) # numero de columnas
    m=length(A[:,1]) # numero de filas
    i=1
    j=0
    while (j<1 && i<=m)
            if (A[i,:]'*x<b[i])
                i=i+1
            else
                j=1
                return A[i,:]
            end
            if (i==m+1)
                return 0
            end
    end
end

function elipsoidal(A,b,D,x,V,v,n)
    i=0
    t=2*(n+1)*(log(V)-log(v))
    while i<t
        i=i+1
        if(pertenece(A,b,x)==0)
            return x
            i=t
        else
            d=pertenece(A,b,x)
            x=eli_vector(x,D,d)
            D=eli_matr(x,D,d)
        end
    end
end


function resolve(A,b,D_0,x_0,x_t,V,v,n,c,c_to_matrix)
    t = 0
    while t<11
        x_t = elipsoidal(A,b,D_0,x_0,V,v,n) # x_t pertenece a la región factible
        b = [b; dot(-c,x_t)] # expandimos el vector b para la nueva restricción
        A = [A; -c_to_matrix]  # expandimos la matriz A para la nueva restricción
        println(x_t)
        println(b)
        println(A)
        println("\n")
        t+=1
    end
    return x_t
end


#Modelo de Dovetail
A = [[1.0, 3.0, 1.0, 0.0] [1.0, 1.0, 0.0, 1.0]]
b = [9.0, 18.0, 7.0, 6.0]
x_0 = [6.0, 6.0] #este punto esta fuera de la región factible
D_0 = [[9.0, 0.0] [0.0, 9.0]] #considerando el radio mayor a 6*sqrt(2) para que envuelva la región factible
V = 343.0 # |x_i| <= 7 => V = 7*7*7 =343
v = (1.0/100.0)
c_to_matrix = [[3.0] [2.0]]
c = [3.0, 2.0]
n = 2.0 #numero de variables
x_t = x_0

println("Loading ...\n")
x_t = resolve(A,b,D_0,x_0,x_t,V,v,n,c,c_to_matrix)
println("Valor objetivo : ", dot(x_t,c))
