function model = changeTFArxnBounds(model,rxnNameList,value,boundType)
	
	%INPUTS
	% model         TFA model structure
	% rxnNameList   List of reactions (cell array or string)
	% value         Bound values
	%               Can either be a vector or a single scalar value if the same
	%               bound value is to be assinged to all reactions
	%
	%OPTIONAL INPUT
	% boundType     'u' - upper, 'l' - lower, 'b' - both 
	%
	%OUTPUT
	% model         TFA model structure with modified reaction bounds
	%
	
	if (nargin < 4)
	    boundType = 'b';
	end
	
	if ~iscell(rxnNameList)
	    rxnNameList = {rxnNameList};
	end
	
	if ~iscell(boundType)
	        boundType = {boundType};
	end
	
	for i = 1:length(rxnNameList) 
	    ind_F = find(ismember(model.varNames,strcat('F_',rxnNameList(i))));
	    ind_R = find(ismember(model.varNames,strcat('R_',rxnNameList(i))));
	
	    switch boundType{i}
	        case 'u'
	            
	            if value(i)>0
	                model.var_lb(ind_F) = 0;
	                model.var_ub(ind_F) = value(i);
	                model.var_lb(ind_R) = 0;
	            elseif value(i)==0
	                model.var_ub(ind_F) = value(i);
	                model.var_lb(ind_F) = value(i);
	                model.var_lb(ind_R) = 0;
	            elseif value(i)<0 && model.var_ub(ind_R) >= abs(value(i))
	                model.var_ub(ind_F) = 0;
	                model.var_lb(ind_F) = 0;
	                model.var_lb(ind_R) = abs(value(i)); 
	            elseif value(i)<0 && model.var_ub(ind_R) < abs(value(i))
	                error('lower bound bigger than upper bound')
	            end
	            
	        case 'l'
	            
	            if value(i)<0
	                model.var_lb(ind_R) = 0;
	                model.var_ub(ind_R) = abs(value(i));
	                model.var_lb(ind_F) = 0;
	            elseif value(i)==0
	                model.var_ub(ind_R) = value(i);
	                model.var_lb(ind_R) = value(i);
	                model.var_lb(ind_F) = 0;
	            elseif value(i)>0 && abs(value(i)) <= model.var_ub(ind_F)
	                model.var_ub(ind_R) = 0;
	                model.var_lb(ind_R) = 0;
	                model.var_lb(ind_F) = value(i);       
	            elseif value(i)>0 && abs(value(i)) > model.var_ub(ind_F)
	                error('lower bound bigger than upper bound')
	            end
	            
	        case 'b'
	            
	           if value(i)>0
	                model.var_ub(ind_F) = value(i);
	                model.var_lb(ind_F) = value(i);
	                model.var_ub(ind_R) = 0; 
	                model.var_lb(ind_R) = 0; 
	            elseif value(i)==0
	                model.var_ub(ind_F) = 0;
	                model.var_lb(ind_F) = 0;
	                model.var_ub(ind_R) = 0; 
	                model.var_lb(ind_R) = 0;
	            elseif value(i)<0
	                model.var_ub(ind_F) = 0; 
	                model.var_lb(ind_F) = 0;  
	                model.var_ub(ind_R) = abs(value(i));
	                model.var_lb(ind_R) = abs(value(i));
	            end
	            
	    end
	end

