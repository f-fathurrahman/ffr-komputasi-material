clear all
clc
tic
% solve the eigenvalue problem
% (£­Delta+V_{ext})u=\lambda u
V=@V_osc_2D;
 %(2) V=@V_r4_2D;
 %(3) V=@V_cos_2D;
 %(4) V=@V_exp2_2D;
 %(5) V=@V_rneg1_2D;
 %(6) V=@V_lnr_2D;
 %(7) V=@V_gauss_2D;
% use the solution (lambda,phi) with Ecx=20 as reference 
% Ec=1:10, plot graphs of residual and |\lambdaN-\lambda| 
Ecx=20;
L=5;
Neig=5;
MAX_Ec=80;
Ecm=10;
%lambdadown_1=2;
%lambdaup_1=4;
%lambdadown_2=4;
%%
[lambda,phi,err_post,dof,vec_k,H,M,vext,HV,u,weight]= eigen_planewave_2D(V,L, Ecx, Neig, MAX_Ec);
%%
lambdaN=zeros(Ecm,Neig);
Delta=zeros(Ecm,Neig);
err=zeros(Ecm,Neig);
erro=zeros(Ecm,Neig);
for Ec=1:Ecm
    [q,p,r,dof0,vec0_k,H0,M0,vext0,HV0,u0,weight0]= eigen_planewave_2D(V,L, Ec, Neig, MAX_Ec);
    %u=real(ifft(phi_m0));
    %[x0,y0] = meshgrid(-5:5/c:5);
    %plot3(x0,y0,u);
    %hold on
    %b=trapz(x,v);
    %I0 = trapz(-5:5/c:5,trapz(-5:5/c:5,u,2));
    %if I0>0
    %p=p;
    %else
    %p=-p;
    %end
    for ll=1:Neig
    lambdaN(Ec,ll)=q(ll);
    %erro(Ec,ll)=trapz(-5:5/(2*MAX_Ec):5,trapz(-5:5/(2*MAX_Ec):5,abs(abs(u0(:,:,ll)).^2-abs(u(:,:,ll)).^2),2));
    err(Ec,ll)=r(ll);
    Delta(Ec,ll)=lambdaN(Ec,ll)-lambda(ll);
    %upper residual
    %cs_ih=lambdadown_2^(-1/2).*(1-lambdaN(:,ll)./lambdadown_2).^(-1);
    %aphih=sqrt(2)*cs_ih.*sqrt(err(:,ll));
    %sinta=err(:,ll)+2*lambdaN(:,ll).*aphih.^2;
    %lower
    %dih=lambdadown_1^2+4*lambdadown_1.*err(:,ll);
    %sonta=max((-lambdaN(:,ll).*aphih.^2+0.5.*(sqrt(dih)-lambdaN(:,ll))),0);
    end
end
toc
%%
ll=1;
semilogy([1:Ecm],Delta(:,ll),'r*--');
axis([1,Ecm,10e-10,1]);
hold on 
semilogy([1:Ecm],err(:,ll),'b*-');
axis([1,Ecm,10e-10,1]);
%hold on
%semilogy([1:Ecm],erro(:,ll).^2,'y*-');
%axis([1,Ecm,10e-8,1]);
%hold on
%semilogy([1:10],sinta,'*-');
%axis([1,10,10e-8,1]);
%hold on
%semilogy([1:10],sonta,'*-');
%axis([1,10,10e-8,1]);
%hold on
legend('|\lambda_{N}-\lambda|','Res');