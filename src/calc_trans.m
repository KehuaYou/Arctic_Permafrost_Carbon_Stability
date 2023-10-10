function calc_trans(i,j)
global COL
global cnstx  cnstTx
global twx twx1 tgx tgx1  ttx
global krw krg
global vsw dnw
global vsg dng
global lambda
global g
global pw pcgw dpth

if i < COL
    krup=2*krw(i,j)*krw(i+1,j)/(krw(i,j)+krw(i+1,j));
    visav=(vsw(i+1,j)+vsw(i,j))/2.0;
    denav=(dnw(i+1,j)+dnw(i,j))/2.0;
    twx(i,j)=cnstx(i,j)*krup/visav;  
    twx1(i,j)=twx(i,j)*denav*g;  


    if pw(i+1,j)+pcgw(i+1,j) > pw(i,j)+pcgw(i,j)+(dng(i+1,j)+dng(i,j))/2*(dpth(i+1,j)-dpth(i,j))*g
        krup=krg(i+1,j);
    else
        krup=krg(i,j);
    end
    visav=(vsg(i+1,j)+vsg(i,j))/2.0;
    denav=(dng(i+1,j)+dng(i,j))/2.0;
    tgx(i,j)=cnstx(i,j)*krup/visav;
    tgx1(i,j)=tgx(i,j)*denav*g;
    lambdaav=(lambda(i+1,j)+lambda(i,j))/2;
    ttx(i,j)=cnstTx(i,j)*lambdaav;
end

if i > 1
    krup=2*krw(i,j)*krw(i-1,j)/(krw(i,j)+krw(i-1,j));
    visav=(vsw(i-1,j)+vsw(i,j))/2.0;
    denav=(dnw(i-1,j)+dnw(i,j))/2.0;
    twx(i-1,j)=cnstx(i-1,j)*krup/visav;
    twx1(i-1,j)=twx(i-1,j)*denav*g;

    if pw(i-1,j)+pcgw(i-1,j) > pw(i,j)+pcgw(i,j)+(dng(i-1,j)+dng(i,j))/2*(dpth(i-1,j)-dpth(i,j))*g
        krup=krg(i-1,j);
    else
        krup=krg(i,j);
    end
    visav=(vsg(i-1,j)+vsg(i,j))/2.0;
    denav=(dng(i-1,j)+dng(i,j))/2.0;
    tgx(i-1,j)=cnstx(i-1,j)*krup/visav;
    tgx1(i-1,j)=tgx(i-1,j)*denav*g;

    lambdaav=(lambda(i-1,j)+lambda(i,j))/2;
    ttx(i-1,j)=cnstTx(i-1,j)*lambdaav;
end




