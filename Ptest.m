
%__________________________________________________________________ %
%                          Multi-Objective                          %
%        Crystal Structure Algorithm (CryStAl) (MOCryStAl)          %
%                                                                   %
%                                                                   %
%                  Developed in MATLAB R2021a (MacOs)               %
%                                                                   %
%                      Author and programmer                        %
%                ---------------------------------                  %
%                      Nima Khodadadi (ʘ‿ʘ)                         %
%                       Siamak Talatahari                           %
%                         Mahdi Azizi                               %
%                         Pooya Sareh                               %
%                                                                   %
%                             e-Mail                                %
%                ---------------------------------                  %
%                         inimakhan@me.com                          % 
%                                                                   %
%                            Homepage                               %
%                ---------------------------------                  %
%                    https://nimakhodadadi.com                      %
%                                                                   %
%                                                                   %
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ----------------------------------------------------------------------- %

function fobj = Ptest(name)

switch name
   
    case 'P8'
        fobj = @P8;
   
end
end

 
% 
% 
%   %% SRN-P8
% % x and y are columnwise, the imput x must be inside the search space and
% % it could be a matrix
function o=P8(x)
o = [0, 0];
%     %SRN
o(1)=(x(1)-2)^2+(x(2)-1)^2+2;
o(2)=9*x(1)-(x(2)-1)^2;
o=o+getnonlinear(x);

    function Z=getnonlinear(x)
        Z=0;
        % Penalty constant
        lam=10^15;
        %        %SRN
        g(1)=x(1)^2+x(2)^2-225;
        g(2)=x(1)-(3*x(2))+10;
        % No equality constraint in this problem, so empty;
        geq=[];
        
        % Apply inequality constraints
        for k=1:length(g)
            Z=Z+ lam*g(k)^2*getH(g(k));
        end
        % Apply equality constraints
        for k=1:length(geq)
            Z=Z+lam*geq(k)^2*getHeq(geq(k));
        end
    end
% Test if inequalities hold
% Index function H(g) for inequalities
    function H=getH(g)
        if g<=0
            H=0;
        else
            H=1;
        end
    end
% Index function for equalities
    function H=getHeq(geq)
        if geq==0
            H=0;
        else
            H=1;
        end
    end
end
%     % ----------------- end ------------------------------
