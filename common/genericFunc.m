function f=genericFunc(d,sfunc)
  sf=['@(x1',sprintf(',x%d',2:d),') ',sfunc,';'];
  f=eval(sf);
end