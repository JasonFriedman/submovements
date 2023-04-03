% cleanExpression - this is a helper function for derivegradianhessian

function cleaned = cleanExpression(s)

thestring = func2str(matlabFunction(s));
closebrackets = strfind(thestring,')');
cleaned = thestring(closebrackets(1)+1:end);