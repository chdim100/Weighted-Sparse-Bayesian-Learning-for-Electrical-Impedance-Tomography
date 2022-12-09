function xin=dual_variable_scaling_rule(xplusdx)
fun=norm(xplusdx);
l=1/fun;
xin=l*xplusdx;
end