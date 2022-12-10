import numpy as np 
import numpy.matlib as mlib
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
#from mpl_toolkits.mplot3d import axes3d
import scipy.io as sio
from scipy.sparse import diags
from scipy.linalg import toeplitz

#%matplotlib qt5


def NOSER_EIT_Reconstruction(J,dv,Lambda):
    ###### J'*J Hessian matrix
    JtJ=J.T.dot(J)
    ###### discrete NOSER filter matrix
    Q=np.diag(np.diag(JtJ))
    dsigma=np.linalg.lstsq(JtJ + Lambda**2 *Q, J.T.dot(dv),rcond=None)
    return dsigma[0]

def EIT_Plot_Image(xc,yc,dsigma,titletext):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    rainb=plt.colormaps['rainbow']
    ax.scatter(xc, yc, dsigma,c = dsigma, marker='o',s=28,cmap=rainb)
    ax.view_init(15, 85)
    ax.set_xlabel('X (A.U.)')
    ax.set_ylabel('Y (A.U.)')
    ax.set_zlabel('$\delta\sigma$ (S/m)')
    plt.title(titletext)
    colnorm = Normalize(vmin=min(dsigma), vmax=max(dsigma)) 
    cbar=plt.colorbar(ScalarMappable(norm=colnorm, cmap=rainb), ax=ax)
    cbar.set_label('$\delta\sigma$ (S/m)', rotation=270,size=13)
    plt.show()
    return(1)


def Estimate_Weights(xc,yc,g,h,naive_ds,std):
    dsr=mlib.repmat(np.abs(naive_ds),1,h)
    wbar=np.zeros((g,1))
    for i in range(g):
        ####### xi belonging in the i-th cluster
        xi,yi=xc[i:i+h],yc[i:i+h]
        di=Eucleidean_dist_square(xc,yc,xi,yi,h)
        Gi=get_Gaussian(di,std,dsr)
        wbar[i]=np.sum(Gi)
    return(wbar/max(wbar))
    
    
def Eucleidean_dist_square(xc,yc,xi,yi,h):
    Lx=np.size(xc)
    xcreformed=np.transpose(mlib.repmat(xc,h,1))
    ycreformed=np.transpose(mlib.repmat(yc,h,1))
    xireformed=mlib.repmat(xi, Lx, 1)
    yireformed=mlib.repmat(yi, Lx, 1)
    return np.power(xcreformed-xireformed,2)+np.power(ycreformed-yireformed,2)

def get_Gaussian(di,std,dsr):
    return(np.multiply(dsr,np.exp(-di/(2*std**2))))
    
    


def Weight_Bound_Opt_SBL(xc,yc,J,dv,h,thetamax,dsigma,std):
    print('Weighted Sparse Bayesian Learning (SBL-weighted):\n')
    M,L=np.shape(J)
    ######Total number of blocks-clusters (aliasing)
    g=L-h+1
    ####Form sparse matrix Psi
    Psi=diags(np.ones(h),0,shape=(L,h)).toarray()
    for i in range(1,g): 
        Psi=np.concatenate((Psi,diags(np.ones(h),-i,shape=(L,h)).toarray()),axis=1)
        #####TO DO: ZERO Initialization of Psi, then complete it block-by-block
    #### Form Phi
    Phi=J.dot(Psi)
    #### Initialization
    theta=1;
    mu=np.zeros((g*h,1))
    gi=np.ones((g,1))
    thresholdgamma=5*1e-02
    w=Estimate_Weights(xc,yc,g,h,dsigma,std)
    Gamma0=0.01*np.sum(np.sqrt(1/(M-1)*np.sum(np.abs(np.power(dv-np.mean(dv),2)),axis=0)))
    rtop=np.zeros((h))
    B=np.zeros((g,h,h))
    So=np.zeros((g*h,g*h))
    Sinv=np.zeros((g*h,g*h))
    for j in range(h):
        rtop[j]=0.9**j
    for i in range(g):
        B[i]=toeplitz(rtop)
        So[h*i:(i+1)*h,:]=np.concatenate((np.zeros((h,i*h)), gi[i]*B[i],np.zeros((h,h*(g-1)-i*h))),axis=1)
        Sinv[h*i:(i+1)*h,:]=np.concatenate((np.zeros((h,i*h)), 1/gi[i]*np.linalg.inv(B[i]),np.zeros((h,h*(g-1)-i*h))),axis=1)
    Btilde=B
    e,emin=1,1e-02
    ####### LOOP
    while e>emin and theta<=thetamax:
        SoPhiT=So.dot(Phi.T)
        PhiSoPhiT=Phi.dot(SoPhiT)
        Sv=Gamma0*np.eye(M)+PhiSoPhiT
        Svinv=np.linalg.inv(Sv)
        muprev=mu
        mu=(SoPhiT.dot(Svinv)).dot(dv)
        Sx=So-(SoPhiT.dot(Svinv)).dot(Phi.dot(So))
        summary=0
        giprev=gi
        for i in range(g):
            Phi_i=Phi[:,i*h:(i+1)*h]
            Sxi=Sx[i*h:(i+1)*h,i*h:(i+1)*h]
            muxi=mu[i*h:(i+1)*h]
            if giprev[i]>thresholdgamma:
                Ai=muxi.T.dot(np.linalg.inv(B[i])).dot(muxi)
                Bhi=np.trace(Svinv.dot(Phi_i).dot(B[i]).dot(Phi_i.T))
                ai=1/w[i]
                bi=w[i]
                zi=(2*ai+np.sqrt(4*ai**2+4*Bhi*(Ai+2*bi)))/(2*(Ai+2*bi))
                gi[i]=1/zi
            else:
                gi[i]=giprev[i]
            summary+=np.trace(Sxi.dot(Phi_i.T).dot(Phi_i))
            Btilde[i]+=1/gi[i]*(Sxi+muxi.dot(muxi.T))
            rtildei=np.mean(np.diag(Btilde[i],1))/np.mean(np.diag(Btilde[i]))
            ri=np.sign(rtildei)*min([abs(rtildei), 0.99])
            rtop2=np.zeros((h))
            for j in range(h):
                rtop2[j]=ri**j
            B[i]=toeplitz(rtop2)
            So[h*i:(i+1)*h,:]=np.concatenate((np.zeros((h,i*h)), gi[i]*B[i],np.zeros((h,h*(g-1)-i*h))),axis=1)
            Sinv[h*i:(i+1)*h,:]=np.concatenate((np.zeros((h,i*h)), 1/gi[i]*np.linalg.inv(B[i]),np.zeros((h,h*(g-1)-i*h))),axis=1)
        ####Update Go 
        Gamma0=1/M*(summary+np.linalg.norm(dv-Phi.dot(mu))**2)
        e=np.linalg.norm(mu-muprev)/np.linalg.norm(mu)
        print(f'Iteration {theta:2.0f} gives error of {e:2.4f}\n')
        theta+=1
    return (Psi.dot(mu))
        
#######load demo model (Adjacent current and voltage pattern, circular model, current
####### amplitude of 1mA p-p)

#####load inverse model pixel grid (circular demo)
xc,yc = np.squeeze(sio.loadmat('demo_inverse_model/coordinates.mat')['xc']),\
    np.squeeze(sio.loadmat('demo_inverse_model/coordinates.mat')['yc'])
    
#####load Jacobian linearized system matrix
J=sio.loadmat('demo_inverse_model/Jacobi.mat')['Jpixel']

##### load demo difference EIT measurements (SNR=40dB)
dv=sio.loadmat('demo_inverse_model/Dv.mat')['dv']

##### NOSER Regularization hyperparameter
Lambda=2.5
#### get naive estimation of conductvity change
ds=NOSER_EIT_Reconstruction(J,dv,Lambda)
#### perform WBO-BSBL
h,thetamax,std=[4,6,(yc[1]-yc[0])/2]
####
deltasigma=Weight_Bound_Opt_SBL(xc,yc,J,dv,h,thetamax,ds,std)
####
EIT_Plot_Image(xc,yc,ds,'NOSER Single-Step')
EIT_Plot_Image(xc,yc,deltasigma,'Weighted Bound Optimization SBL')
