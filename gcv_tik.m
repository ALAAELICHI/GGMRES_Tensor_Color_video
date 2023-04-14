function alpha = gcv_tik(f,g,betta)
n=size(f,1);
t=f(n);
m=size(g,1);
f=f(1:m,:);
alpha = fminbnd(@(alpha) GCV(alpha,f,g),0,0.0001);

    function G = GCV(alpha, f, g)
        


u=((alpha^2*f)./(g.^2 + alpha^2)).^2;

v=alpha^2./(g.^2 + alpha^2);

G=betta^2*(sum(u)+t^2)/((sum(v)+1)^2);

  end
  
end