function y = ProjectionLS( Pt,b,T, transp_flag )
It = speye(size(Pt,1));

if strcmp(transp_flag,'notransp')   
It = It(:,T);
y = It*b - Pt*(Pt(T,:)'*b);

elseif strcmp(transp_flag,'transp')

It = It(T,:);
y = It*b - Pt(T,:)*(Pt'*b);

end
end

