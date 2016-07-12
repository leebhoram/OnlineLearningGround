% by Bhoram Lee
% Apr 2014
classdef KF1D
    
    properties
        x
        z
        px
        q
        r   
        r_default
    end
    
    methods
        function fp = initialize(fp, x0, q, r)            
            fp.x = x0;
            fp.q = q;
            fp.r = r;
            fp.r_default = r;
            fp.px = 2*q;           
        end
        
        function [fp, x_post, xp_post] = update(fp, meas, r_)
            x_post = 0;
            xp_post = 0;
            if ~isempty(fp.q) && ~isempty(fp.r)
                fp.z = meas;   
                fp = propagate(fp);
                if nargin == 3
                    fp.r = r_; % adaptive noise model
                else
                    fp.r = fp.r_default;
                end
                fp = measUpdate(fp);
                
                x_post = fp.x;
                xp_post = fp.px;
            end
        end
        
        function [fp,x_pred,xp_pred]= propagate(fp)
            % model 
            fp.x = fp.x;
            fp.px = fp.px + fp.q;
            x_pred = fp.x;
            xp_pred = fp.px;
        end
        
         function [fp, x_post, xp_post] = zero(fp)
            x_post = 0;
            xp_post = 0;
            if ~isempty(fp.q) && ~isempty(fp.r)
                fp.x = 0;
                fp.px = fp.q;         
            end
        end
    end
    
    methods (Access = private)
                
        function fp = measUpdate(fp)    
            
            sum_sx_sy = (fp.px + fp.r);
            fp.x = (fp.r*fp.x + fp.px*fp.z)/sum_sx_sy;
            fp.px = fp.r*fp.px/sum_sx_sy;                 
        end
    end
    
end

