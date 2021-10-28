function [resstp1]  =   solverstp1(Hez,sz,gam_try,t)
    % "yalmip('clear')" is a good practive for iterative purposes.
    % Although new variables are assigned at each iteration, old variables
    % grows the stored variables and it affects solver somehow. 
    yalmip('clear');
    g   = gam_try;
    
    nw  =   sz.nw;
    ne  =   sz.ne;
    na  =   sz.na;
    nx  =   sz.nx;
    
    % Lyap supply var and supply func for interconnections (na many dec var matrices)
    zzz11 = [];
    zzz22 = [];
    zzz12 = [];
    xx  = [];
    for i = 1:na % it was na before
        % X is the decision variable for Lyap. Supply Func.
        % X((i-1)*nx+1:i*nx,(i-1)*nx+1:i*nx) = sdpvar(nx, nx, 'symmetric');
        x           =   sdpvar(nx, nx, 'symmetric');
        x           =   reshape(x,1,[]);
        xx          =   [xx,x];
        % X11 and X12 are decision variables for neutral interconnections
%         X11((i-1)*nw+1:i*nw,1:nw*2) =   [sdpvar(nw,nw,'symmetric'),sdpvar(nw,nw,'symmetric')];
        X11((i-1)*nw+1:i*nw,1:nw*2) =   [sdpvar*eye(nw),sdpvar*eye(nw)];
        x12                         =   sdpvar(nw,nw,'skew','real');
        X12((i-1)*nw+1:i*nw,1:nw*2) =   [x12,x12'];
        z11                         =   blkdiag(zeros(na*nw),-X11((i-1)*nw+1:i*nw,1:nw));
        z22                         =   blkdiag(zeros(na*nw), X11((i-1)*nw+1:i*nw,nw+1:2*nw));
        z12                         =   blkdiag(zeros(na*nw), X12((i-1)*nw+1:i*nw,nw+1:2*nw));

        % Vectorize z11, z12 and z22 to be able to append proceeding variables
        % for all agents
        % This is required due to how "sdpvar" is defined in YALMIP.

        zz11    =   reshape(z11,1,((na+1)*nw)^2);
        zz22    =   reshape(z22,1,((na+1)*nw)^2);
        zz12    =   reshape(z12,1,((na+1)*nw)^2);

        zzz11     =   [zzz11,zz11];
        zzz22     =   [zzz22,zz22];
        zzz12     =   [zzz12,zz12];
        %
        %

        % Final supply function for virtual agent that contains interconnection
        % info is written in the following script
        
        if i == na
            z11 = []; z12 = []; z22 = [];
            for ii = 1:na
                z11 = blkdiag(z11,-X11((ii-1)*nw+1:ii*nw,nw+1:2*nw));
                z22 = blkdiag(z22, X11((ii-1)*nw+1:ii*nw,1:nw));
                z12 = blkdiag(z12,-X12((ii-1)*nw+1:ii*nw,1:nw));
                if ii == na
                    z11 = blkdiag(z11,zeros(nw));
                    z22 = blkdiag(z22,zeros(nw));
                    z12 = blkdiag(z12,zeros(nw));
                end
            end
            zz11    =   reshape(z11,1,((na+1)*nw)^2);
            zz22    =   reshape(z22,1,((na+1)*nw)^2);
            zz12    =   reshape(z12,1,((na+1)*nw)^2);

            zzz11     =   [zzz11,zz11];
            zzz22     =   [zzz22,zz22];
            zzz12     =   [zzz12,zz12];
        end
    end
    X   =   reshape(xx,nx,nx,[]);
    Z11 =   reshape(zzz11,nw*ne,nw*ne,[]);
    Z22 =   reshape(zzz22,nw*ne,nw*ne,[]);
    Z12 =   reshape(zzz12,nw*ne,nw*ne,[]);
    Z21 =   Z12';
    % In the following for loop, LMI is converted to full rank.
    for i = 1:ne
        if i<ne
        C   =   [X(:,:,i)*Hez(i).A+Hez(i).A'*X(:,:,i), X(:,:,i)*Hez(i).B;
                 Hez(i).B'*X(:,:,i), -Hez(i).I] + ...
                 1/(g^2)*[Hez(i).C, Hez(i).D]'*[Hez(i).C, Hez(i).D] + ...
                 [Hez(i).C, Hez(i).D; zeros(nw*ne,nx), eye(nw*ne)]'*[Z11(:,:,i),Z12(:,:,i);Z21(:,:,i),Z22(:,:,i)]*[Hez(i).C, Hez(i).D; zeros(nw*ne,nx), eye(nw*ne)];

        nc      =  size(C,1);
        nzeros  = 0;
        for j=1:nc
            if isempty(C(j,j))
                nzeros = nzeros + 1;
                ind(1,nzeros) = j;
            elseif value(C(j,j))==0
                nzeros = nzeros + 1;
                ind(1,nzeros) = j;
            end
        end
        C(:,ind)    = [];
        C(ind,:)    = [];
        LMI(i).C    = C;
        C   = [];
        ind = [];
        else
            C   =   [Hez(i).A+Hez(i).A', Hez(i).B;
                    Hez(i).B', -Hez(i).I] + ...
                     1/(g^2)*[Hez(i).C, Hez(i).D]'*[Hez(i).C, Hez(i).D] + ...
                     [Hez(i).C, Hez(i).D; zeros(nw*ne,nx), eye(nw*ne)]'*[Z11(:,:,i),Z12(:,:,i);Z21(:,:,i),Z22(:,:,i)]*[Hez(i).C, Hez(i).D; zeros(nw*ne,nx), eye(nw*ne)];

            nc      =  size(C,1);
            nzeros  = 0;
            for j=1:nc
                if isempty(C(j,j))
                    nzeros = nzeros + 1;
                    ind(1,nzeros)   = j;
                elseif value(C(j,j))==0
                    nzeros = nzeros + 1;
                    ind(1,nzeros)   = j;
                end
            end
            C(:,ind)    = [];
            C(ind,:)    = [];
            LMI(i).C    = C-t*eye(size(C));
            C       = [];
            ind     = [];
            nzeros  = [];
            nc      = [];
        end
    end

    % Other constraints on the optimization is defined in the following for
    % loop
        for i = 1:ne-1
            LMI(i+ne).C    =   X(:,:,i);
        end
    
    % Collect all the constraints for YALMIP
%     Cons       = [LMI(1).C<=0;LMI(2).C<=0;LMI(3).C<=0;LMI(4).C<=0;LMI(5).C<=0;
%                 LMI(6).C>=0;LMI(7).C>=0;LMI(8).C>=0;LMI(9).C>=0];
    % one can collect constraints iteratively for YALMIP as follows
    sCons   =   size(LMI,2);
    Cons    =   [];
    for i = 1: sCons
        if i<=sCons-na
            Cons    =   [Cons, (LMI(i).C <= 0):['Inequality Constraint  ' num2str(i)]];
        else
            Cons    =   [Cons, (LMI(i).C >= 0):['Lyapunov Storage  ' num2str(i)]];
        end
    end
    opt         =	sdpsettings('verbose',0,'warning',0);
    solTime     =   cputime;
    solndiag    =   optimize(Cons,[],opt);
    solTime     =   cputime - solTime;
    Xval        =   value(X);

    % Some of these matrices Z are NAN as they disappear while the constraints
    % are calculated. Though we still extract their values for generality of
    % the solver.
    if solndiag.problem~=0
        resstp1.succeed = false;
    elseif solndiag.problem==0
        resstp1.succeed = true;
        for i = 1:na
            X11cval(:,:,i)      =   value(X11((i-1)*nw+1:i*nw,nw+1:2*nw));
            X11ival(:,:,i)      =   value(X11((i-1)*nw+1:i*nw,1:nw));
            X12cval(:,:,i)      =	value(X12((i-1)*nw+1:i*nw,nw+1:2*nw));
            lmicheck(:,:,i)     =	value(LMI(i).C);
            max_eigenLMI(:,:,i)	=   max(eig(lmicheck(:,:,i)));
            min_eigenLMI(:,:,i)	=   min(eig(lmicheck(:,:,i)));
            max_eigenX(:,:,i)	=   max(eig(Xval(:,:,i)));
            min_eigenX(:,:,i)	=   min(eig(Xval(:,:,i)));
            EigX(:,:,i)         =	eig(Xval(:,:,i));
            EucNorm(:,:,i)      =   norm(EigX(:,:,i),2);
        end
        if max(max_eigenLMI) < 0 && min(min_eigenLMI) < 0
            resstp1.X       =   Xval;
            resstp1.X11c    =   X11cval;
            resstp1.X11i    =   X11ival;
            resstp1.X12c    =   zeros(size(X12cval));
            resstp1.g       =   g;
            resstp1.Z11     =   [];
            resstp1.Z22     =   [];
            resstp1.Z12     =   [];
            for i = 1:na
                resstp1.Z11(:,:,i)     =   blkdiag(zeros(na*nw),-resstp1.X11i(:,:,i));
                resstp1.Z22(:,:,i)     =   blkdiag(zeros(na*nw), resstp1.X11c(:,:,i));
                resstp1.Z12(:,:,i)     =   blkdiag(zeros(na*nw), resstp1.X12c(:,:,i));
                resstp1.Z(:,:,i)       =   [resstp1.Z11(:,:,i),resstp1.Z12(:,:,i);resstp1.Z12(:,:,i)',resstp1.Z22(:,:,i)];
                if i == na
                    z11 = []; z12 = []; z22 = [];
                    for ii = 1:na
                        z11 = blkdiag(z11,-resstp1.X11c(:,:,ii));
                        z22 = blkdiag(z22, resstp1.X11i(:,:,ii));
                        z12 = blkdiag(z12,-resstp1.X12c(:,:,ii));
                        if ii == na
                            z11 = blkdiag(z11,zeros(nw));
                            z22 = blkdiag(z22,zeros(nw));
                            z12 = blkdiag(z12,zeros(nw));
                        end
                    end
                    resstp1.Z11(:,:,i+1)     =   z11;
                    resstp1.Z22(:,:,i+1)     =   z22;
                    resstp1.Z12(:,:,i+1)     =   z12;
                    resstp1.Z(:,:,i+1)       =   [resstp1.Z11(:,:,i+1),resstp1.Z12(:,:,i+1);resstp1.Z12(:,:,i+1)',resstp1.Z22(:,:,i+1)];
                end
            end
            resstp1.Hez =   Hez;
            resstp1.Adj =   Hez(1,ne).D(na*nw+1:ne*nw,na*nw+1:ne*nw);
            % resstp1.Z(:,:,5) which represents the final inequalities Z matrix is
            % denoted as Zc. Then Zc is described with a definition as 
            % Zc = vc_mod*e_mod*vc_mod'
%             Zc      = resstp1.Z(:,:,5);
%             [vc,ec] = eig(Zc);
%             ecv     = diag(ec);
%             vc_mod = vc*diag(sqrt(abs(diag(ec))));
%             for i = 1:length(ecv)
%                 if (ecv(i,1))~=0
%                     zcv(i,1)=sign(ecv(i,1));
%                 else
%                     zcv(i,1)=0;
%                 end
%             end
%             e_mod   =   diag(zcv);
%             resstp1.Zc = Zc;
%             resstp1.vc = vc_mod;
%             resstp1.ec = e_mod;
        end
    end
end