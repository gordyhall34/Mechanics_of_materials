import math
import matplotlib.pyplot as plt
import numpy as geek


def main():
    sig_x, sig_y, sig_z, tau_xy, tau_yz, tau_zx = values()
    calculations(sig_x, sig_y, sig_z, tau_xy, tau_yz, tau_zx)


def values():
    sig_x1 = input("Please input the Sigma X value: ")
    sig_x = float(sig_x1)
    print(sig_x)

    sig_y1 = input("Please input the Sigma Y value: ")
    sig_y = float(sig_y1)
    print(sig_y)
    sig_z1 = input("Please input the Sigma Z value: ")
    sig_z = float(sig_z1)
    print(sig_z)
    tau_xy1 = input("Please input the Tau XY value: ")
    tau_xy = float(tau_xy1)
    print(tau_xy)
    tau_yz1 = input("Please input the Tau YZ value: ")
    tau_yz = float(tau_yz1)
    print(tau_yz)
    tau_zx1 = input("Please input the Tau ZX value: ")
    tau_zx = float(tau_zx1)
    print(tau_zx)
    return sig_x, sig_y, sig_z, tau_xy, tau_yz, tau_zx


def calculations(sig_x, sig_y, sig_z, tau_xy, tau_yz, tau_zx):
    sig_h = (1 / math.sqrt(2)) * math.sqrt((((sig_x - sig_y) ** 2) + (sig_y - sig_z) ** 2) + (sig_z - sig_x) ** 2) + (
            6 * ((tau_xy ** 2) + (tau_yz ** 2) + (tau_zx ** 2)))
    print(sig_h)

    ## Principle Strains
    i_1 = sig_x + sig_y + sig_z

    i_2 = (sig_x * sig_y) + (sig_y * sig_z) + (sig_z * sig_x) - (tau_xy ** 2) - (tau_yz ** 2) - (tau_zx ** 2)

    i_3 = (sig_x * sig_y * sig_z) + (2 * tau_xy * tau_yz * tau_zx) - (sig_x * (tau_yz ** 2)) - (
            sig_y * (tau_zx ** 2)) - (sig_z * (tau_xy ** 2))

    p = [1, -i_1, i_2, -i_3]
    cubic = geek.roots(p)

    b = sorted(cubic, reverse=True)

    sig_1 = b[0]
    print(sig_1)
    sig_2 = b[1]
    print(sig_2)
    sig_3 = b[2]
    print(sig_3)

    ## Principle Stains
    tau_1 = abs(sig_2 - sig_3) / 2  ##Finds Principle Strain 1
    print(tau_1)
    tau_2 = abs(sig_1 - sig_3) / 2  ##Finds Principle Strain 2
    print(tau_2)
    tau_3 = abs(sig_1 - sig_2) / 2  ##Finds Principle Strain 3
    print(tau_3)

    ## MAX Ultimate stress C1
    sig_1abs = abs(sig_1)  ##Makes the principle stresses aboslute values
    sig_2abs = abs(sig_2)  ##Makes the principle stresses aboslute values
    sig_3abs = abs(sig_3)  ##Makes the principle stresses aboslute values
    array = [sig_1abs, sig_2abs, sig_3abs]  ##puts stresses in an array
    sig_n = max(array)  ##Finds max value of array
    sig_u = sig_n * 2  ##Using the given safety factor, finds the ultimate stress

    print(["Using the Maximum Normal Stress Fracture Criterion and a safety factor of 2, the material needs to have "
           "this ultimate tensile strength: ", str(sig_u)])

    ## MAX shear stress C2 need help with this one
    tarray = [tau_1, tau_2, tau_3]  ##makes the principle strains into an array
    e1 = abs(sig_2 - sig_3)  ##Makes the principle strains aboslute values
    e2 = abs(sig_1 - sig_3)  ##Makes the principle strains aboslute values
    e3 = abs(sig_1 - sig_2)  ##Makes the principle strains aboslute values
    sig_sarray = [e1, e2, e3]  ##Places the absolute values of the strains into an array
    tau_yield = max(tarray)  ##Finds the max value of the principle strain array
    sig_s = max(sig_sarray)  ##finds the max value of the absolute value of the principle strain array
    sig_0 = sig_s * 2  ##Using the given safety factor, finds the ultimate
    print(['Using the Maximum Shear Stress Yield Criterion and a safety factor of 2,'
          ' the material needs to have this yield strength : ', str(sig_0)])

    ## OCT shear stress C3
    tau_ho= (1/3)*math.sqrt(((sig_1 - sig_2)**2)+((sig_2 - sig_3)**2)+((sig_1 - sig_3)**2)) ##Finds the octahedral shear strain
    sig_h = (1/math.sqrt(2))*math.sqrt((((sig_x - sig_y)**2)+(sig_y-sig_z)**2)+(sig_z - sig_x)**2)+(6*((tau_xy**2)+(tau_yz**2)+(tau_zx**2))) ##finds the octahedral shear stress
    sig_oct = sig_h*2 ##Multiples the
    print([' Using the Octahedral Shear Stress Yield Criterion and a safety factor of 2, the material needs to have this yield strength: ', str(sig_oct)])

    ## Plotting the 3-D Mohr's Circle
    angles = 0:2*math.pi:0.01
    O_ave1=[(O_1+O_3)/2 0];
    O_ave2=[(O_1+O_2)/2 0];
    O_ave3=[(O_2+O_3)/2 0];
    cirlce1=[O_ave1(1)+t_2*cos(angles') O_ave1(2)+t_2*sin(angles')];
    cirlce2=[O_ave2(1)+t_3*cos(angles') O_ave2(2)+t_3*sin(angles')];
    cirlce3=[O_ave3(1)+t_1*cos(angles') O_ave3(2)+t_1*sin(angles')];
    plot(cirlce1(:,1),cirlce1(:,2),'b',cirlce2(:,1),cirlce2(:,2),'g',...
        cirlce3(:,1),cirlce3(:,2),'r');axis equal;grid on
%{
this code is to build the MOHR's circle on a graph,
specifically the three circles if they are needed
%}
%ANNOTATIONS
    if O_1>0
        text(O_1*1.01,0,'\sigma_1','fontsize',15);
    else
        text(O_1*0.95,0,'\sigma_1','fontsize',15);
    end
    if O_2>0
        text(O_2*1.1,0,'\sigma_2','fontsize',15);
    else
        text(O_2*0.99,0,'\sigma_2','fontsize',15);
    end
    if O_3>0
        text(O_3*1.1,0,'\sigma_3','fontsize',15);
    else
        text(O_3*0.99,0,'\sigma_3','fontsize',15);
    end
    text(O_ave1(1),t_1*0.9,'\tau_1','fontsize',15);
    text(O_ave2(1),t_2*0.9,'\tau_2','fontsize',15);
    text(O_ave3(1),t_3*0.9,'\tau_3','fontsize',15);
    xlabel('Normal Stress, \sigma','fontsize',15);
    ylabel('Shear Stess, \tau','fontsize',15)
%{
This adds Markers to the graph,
so that each circle is labeled
%}
%% E
    S= [O_x  t_xy t_zx;t_xy O_y t_yz;t_zx t_yz O_z];
    [Eigenvectors,Eigenvalues]= eig(S); %finds eigenvectors and egienvalues
    sigmai=Eigenvalues;
    ni= Eigenvectors;
    sigma1=sigmai(1,1);
    sigma2=sigmai(2,2);
    sigma3=sigmai(3,3);
    n1=ni(1,1);
    n2=ni(2,2);
    n3=ni(3,3);
    pair1=[sigma1,n1]
    pair2=[sigma2,n2]
    pair3=[sigma3,n3]
%{
This finds the eigenvalues
and the eigenvecors, places
them into a matrix and
pulls them from their
locations and pairs them
with their matching value
%}
if __name__ == '__main__':
    main()
