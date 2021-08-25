
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
function costs=GetCosts(pop)

    nobj=numel(pop(1).Cost);
    costs=reshape([pop.Cost],nobj,[]);

end