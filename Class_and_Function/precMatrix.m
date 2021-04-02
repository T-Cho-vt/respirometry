classdef precMatrix
    %precMatrix class
    %  A precMatrix object is used to represent a preconditioner matrix P used
    %  for the Respirometry example.
    %
    %  The precMatrix class has inputs:
    %   h               - impulse response function
    %   tau             - tolerance for the preconditioner
    %   and is based on a structure with two fields:
    %       lambda - eigenvalues
    %   transpose - indicates if the matrix has been transposed
    %
    %  Calling Syntax:
    %
    %    P = precMatrix    (returns object with empty fields)
    %    P = precMatrix(precMatrix) (returns precMatrix)
    %    P = precMatrix(h, tau)
    %
    % J. Chung, 5/4/15
    
    properties
        lambda
        transpose
    end
    
    methods
        function P = precMatrix(varargin) % Constructor
            switch nargin
                case 0
                    P.transpose = false;
                    P.lambda = [];
                case 1
                    if isa(varargin{1}, 'precMatrix')
                        P = varargin{1};
                    else
                        error('Incorrect input arguments')
                    end
                otherwise
                    P.transpose = false;
                    if nargin == 2
                        h = varargin{1};
                        tau = varargin{2};
                        lam = fft(h);
                        lam(abs(lam)<tau) = 1;
                        P.lambda = lam;
                    else
                        error('Incorrect input arguments')
                    end
            end
        end
        
        
        function P = ctranspose(P) % Overload transpose
            if P.transpose == 0
                P.transpose = 1;
            else
                P.transpose = 0;
            end
        end
        
        %         function P = mrdivide(P,arg2) % Overload matrix scalar division
        %             P.tau = P.tau/arg2;
        %         end
        
        function y = mldivide(arg1,arg2) % Overload matrix backslash
            %  Check to see if arg2 is not a vector
            if prod(size(arg2)) ~= length(arg2)
                disp('Warning:precMatrix/mtimes assumes input is a vector')
            end
            s = arg2(:);
            lambda = arg1.lambda;
            
            if arg1.transpose
                % Matrix transpose times vector
                y = real(fft(ifft(s)./lambda));
%                 y = (fft(ifft(s)./lambda));
            else
                % Matrix times vector
                y = real(ifft(fft(s)./lambda));
%                 y = (ifft(fft(s)./lambda));
            end
        end
        
        function y = mtimes(arg1, arg2) % Overload matrix vector multiplicaiton
            %   Implement P*s and P'*s for precMatrix object P and "vector" s.
            
            if ( isa( arg1, 'precMatrix' ) )
                % check to see of arg2 is a scalar
                if (( isa ( arg2, 'double' )) && (length(arg2)==1))
                    error('Not implemented yet')
                else
                    %  Check to see if arg2 is not a vector
                    if prod(size(arg2)) ~= length(arg2)
                        disp('Warning:precMatrix/mtimes assumes input is a vector')
                    end
                    s = arg2(:);
                    lambda = arg1.lambda;
                    
                    if arg1.transpose
                        % Matrix transpose times vector
                        y = real(fft(ifft(s).*lambda));
%                         y = (fft(ifft(s).*lambda));
                    else
                        % Matrix times vector
                        y = real(ifft(fft(s).*lambda));
%                         y = (ifft(fft(s).*lambda));
                    end
                end
                
            elseif (( isa(arg1, 'double')) && (length(arg1)==1))
                % check to see of arg1 is a scalar
                %                 y = arg2;
                %                 y.tau = y.tau/arg1;
                
                error('Multiplication is not implemented.')
            else
                error('Multiplication is not implemented.')
            end
        end
        
        function varargout = size( P, dim ) % Overload size for precMatrix object
            d(1) = length(P.h);
            d(2) = length(P.h);
            
            if nargout == 1 || nargout == 0
                if nargin > 1
                    varargout{1} = d(dim);
                else
                    varargout{1} = d;
                end
            else
                varargout{1} = d(1);
                varargout{2} = d(2);
            end
        end
        
    end % methods
end % classdef

