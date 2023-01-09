MCP_func = function(z,lambda,a=3){
  MCP.value = (lambda*abs(z) - z^2/2/a) * (abs(z) <= lambda*a) +
    (a*lambda^2/2) * (abs(z) > lambda*a)
  return(MCP.value)
}
