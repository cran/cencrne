MCP_soft = function(z,lambda,a=3,kappa=1){

  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: MCP_soft
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description:
  ##            Define the MCP threshold operator.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ z: a float value or a vector, the independent variable in the MCP.
  ## @ lambda: a float value, the tuning parameter in the MCP.
  ## @ a: a float value, regularization parameter in the MCP, the default setting is 3.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ The result of the MCP threshold operator.
  ## -----------------------------------------------------------------------------------------------------------------
  norm.z = sqrt(sum(z^2))
  return( S_soft(z,lambda/kappa)/(1-1/a/kappa) * (norm.z - a*lambda <= 0) + z * (norm.z - a*lambda > 0) )
}
