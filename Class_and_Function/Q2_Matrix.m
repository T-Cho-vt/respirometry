classdef Q2_Matrix
    
    properties
        n
        transpose
    end
    
    methods
        function P = Q2_Matrix(varargin) % constructur
            switch nargin
                case 1
                    P.transpose = false;
                    P.n = varargin{1};
                otherwise
                    error('Incorrect number of input arguments')
            end % switch
        end 
        
        function P = ctranspose(P) % Overload transpose
            P.transpose = not(P.transpose); % switches boolean trnaspose flag
        end % transpose
            
        
        function y = mtimes(H,x)
            
            if length(x) ~= H.n
                error('Length of vector must match the number of columns of matrix')
            end
            
            if H.transpose
               y = x- [2*x(2:end); 0] + [x(3:end); 0; 0;]; 
            else
               y = x - [0; 2*x(1:end-1)] +[0; 0; x(1:end-2)]; 
            end
            
        end
        
        function varargout = size(H,dim)
            d(1) = H.n;
            d(2) = H.n;

            if nargout == 1 || nargout == 0
                if nargin >1
                    varargout{1} = d(dim);
                else
                    varargout{1} = d;
                end
            else
                varargout{1} = d(1);
                varargout{2} = d(2);
            end
        end % size
        
        function l = length(H)
            l = H.n;
        end % length
        
    end % methods
    
    
end % classdef