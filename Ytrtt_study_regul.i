/* ÉTUDE DES RÉGULARISATIONS SPATIO-TEMPORELLES
 *
 * Si l'on suppose que la dynamique de notre objet spatio-temporel consiste à
 * voir des fronts se déplacer, alors :
 *
 * -> 1 condition nécessaire : un gradient temporel est forcément lié à un
 * gradient spatial.
 *
 * -> mais pas suffisante : un gradient spatial n'est pas forcément lié à un
 *    gradient temporel.
 *
 * Il faut donc une régularisation qui favorise les gradients temporels aux
 * lieux des gradients spatiaux.
 *
 * On prend le cas simple suivant d'un objet 2D+t :
 *
 * -------------------     -------------------     -------------------
 * |     |     |     |     |     |     |     |     |     |     |     |
 * |  0  |  0  |  0  |     |  0  |  0  |  0  |     |  0  |  0  |  0  |
 * |     |     |     |     |     |     |     |     |     |     |     | 
 * -------------------     -------------------     -------------------
 * |     |     |     |     |     |     |     |     |     |     |     |
 * |  0  |  0  |  0  |     |  ?  |  ?  |  0  |     |  1  |  0  |  0  |
 * |     |     |     |     |     |     |     |     |     |     |     | 
 * -------------------     -------------------     -------------------
 * |     |     |     |     |     |     |     |     |     |     |     |
 * |  0  |  0  |  0  |     |  ?  |  ?  |  0  |     |  1  |  1  |  0  |
 * |     |     |     |     |     |     |     |     |     |     |     | 
 * -------------------     -------------------     -------------------
 *
 *       frame 1                 frame 2                 frame 3
 *
 * Un front se déplace en diagonal à partir du coin inférieur gauche. On cherche
 * à placer un pixel dans le frame 2 pour assurer la continuité temporelle entre
 * les frames 1 et 3. On voit bien qu'il doit idéalement se placer dans la case
 * en bas à gauche pour satisfaire cette condition. Il s'agit donc de trouver la
 * régularisation pour laquelle le critère est le moins coûteux lorsque le pixel
 * est placé à cette position idéale. On teste plusieurs régularisations de type
 * TV :
 *
 * -> TV non isotrope non séparable.
 * -> TV non isotrope séparable.
 * -> TV isotrope non séparable.
 * -> TV isotrope séparable.
 *
 * Notre petite expérience montre que quel que soit le critère utilisé, la
 * position la "moins chère" est bien celle escomptée.
 *
 */

func TV_non_isotrope_non_separable(x)
{
    n=2;
    fx=0;
    for(k=1;k<=n;++k) {
        for(i=1;i<=n;++i) {
            for(j=1;j<=n;++j) {
                gs = abs(x(i,j,k)-x(i+1,j,k))^2 + abs(x(i,j,k)-x(i,j+1,k))^2;
                gt = abs(x(i,j,k)-x(i,j,k+1))^2;
                fx += sqrt(gs+gt);
            }
        }
    }

    return fx;
}

func TV_non_isotrope_non_separable_reverse(x)
{
    n=2;
    fx=0;
    for(k=2;k<=n+1;++k) {
        for(i=1;i<=n;++i) {
            for(j=1;j<=n;++j) {
                gs = abs(x(i,j,k)-x(i+1,j,k))^2 + abs(x(i,j,k)-x(i,j+1,k))^2;
                gt = abs(x(i,j,k-1)-x(i,j,k))^2;
                fx += sqrt(gs+gt);
            }
        }
    }

    return fx;
}


func TV_non_isotrope_separable(x)
{
    n=2;
    fx=0;
    for(k=1;k<=n;++k) {
        for(i=1;i<=n;++i) {
            for(j=1;j<=n;++j) {
                gs = abs(x(i,j,k)-x(i+1,j,k))^2 + abs(x(i,j,k)-x(i,j+1,k))^2;
                gt = abs(x(i,j,k)-x(i,j,k+1));
                fx += sqrt(gs)+gt;
            }
        }
    }

    return fx;
}


func TV_isotrope_non_separable(x)
{
    n=2;
    fx=0;
    for(k=1;k<=n;++k) {
        for(i=1;i<=n;++i) {
            for(j=1;j<=n;++j) {
                gsx = 0.25 * (abs(x(i,j,k)-x(i+1,j,k))^2 + abs(x(i,j+1,k)-x(i+1,j+1,k))^2 + abs(x(i,j,k+1)-x(i+1,j,k+1))^2 + abs(x(i,j+1,k+1)-x(i+1,j+1,k+1))^2);
                gsy = 0.25 * (abs(x(i,j,k)-x(i,j+1,k))^2 + abs(x(i+1,j,k)-x(i+1,j+1,k))^2 + abs(x(i,j,k+1)-x(i,j+1,k+1))^2 + abs(x(i+1,j,k+1)-x(i+1,j+1,k+1))^2);
                gt = 0.25*(abs(x(i,j,k)-x(i,j,k+1))^2 + abs(x(i+1,j,k)-x(i+1,j,k+1))^2 + abs(x(i,j+1,k)-x(i,j+1,k+1))^2 + abs(x(i+1,j+1,k)-x(i+1,j+1,k+1))^2);
                fx += sqrt(gsx + gsy + gt);
            }
        }
    }

    return fx;
}

func TV_isotrope_separable(x)
{
    n=2;
    fx=0;
    for(k=1;k<=n;++k) {
        for(i=1;i<=n;++i) {
            for(j=1;j<=n;++j) {
                gsx = 0.25 * (abs(x(i,j,k)-x(i+1,j,k))^2 + abs(x(i,j+1,k)-x(i+1,j+1,k))^2 + abs(x(i,j,k+1)-x(i+1,j,k+1))^2 + abs(x(i,j+1,k+1)-x(i+1,j+1,k+1))^2);
                gsy = 0.25 * (abs(x(i,j,k)-x(i,j+1,k))^2 + abs(x(i+1,j,k)-x(i+1,j+1,k))^2 + abs(x(i,j,k+1)-x(i,j+1,k+1))^2 + abs(x(i+1,j,k+1)-x(i+1,j+1,k+1))^2);
                gt = 0.25*(abs(x(i,j,k)-x(i,j,k+1))^2 + abs(x(i+1,j,k)-x(i+1,j,k+1))^2 + abs(x(i,j+1,k)-x(i,j+1,k+1))^2 + abs(x(i+1,j+1,k)-x(i+1,j+1,k+1))^2);
                fx += sqrt(gsx + gsy) + sqrt(gt);
            }
        }
    }

    return fx;
}

n=3;
x=array(double,n,n,n);
x(1,1,1)=1.0;
x(1,1,3)=1.0;
x(1,3,3)=1.0;
x(3,1,3)=1.0;
x(2,2,3)=1.0;
x(1,2,3)=1.0;
x(2,1,3)=1.0;

ncase=7;
TV1=array(double,ncase);
TV2=array(double,ncase);
TV3=array(double,ncase);
TV4=array(double,ncase);

// pixel en 1, 2 et 3 => le cas que l'on veut !
x1=x;
x1(1,2,2)=1.0;
x1(2,1,2)=1.0;
x1(1,1,2)=1.0;
TV1(1)=TV_non_isotrope_non_separable(x1);
TV2(1)=TV_non_isotrope_separable(x1);
TV3(1)=TV_isotrope_non_separable(x1);
TV4(1)=TV_isotrope_separable(x1);

// pixel en 1
x2=x;
x2(1,1,2)=1.0;
TV1(2)=TV_non_isotrope_non_separable(x2);
TV2(2)=TV_non_isotrope_separable(x2);
TV3(2)=TV_isotrope_non_separable(x2);
TV4(2)=TV_isotrope_separable(x2);

// pixel en 2
x3=x;
x3(2,1,2)=1.0;
TV1(3)=TV_non_isotrope_non_separable(x3);
TV2(3)=TV_non_isotrope_separable(x3);
TV3(3)=TV_isotrope_non_separable(x3);
TV4(3)=TV_isotrope_separable(x3);

// pixel en 3
x4=x;
x4(1,2,2)=1.0;
TV1(4)=TV_non_isotrope_non_separable(x4);
TV2(4)=TV_non_isotrope_separable(x4);
TV3(4)=TV_isotrope_non_separable(x4);
TV4(4)=TV_isotrope_separable(x4);

// aucun
x5=x;
TV1(5)=TV_non_isotrope_non_separable(x5);
TV2(5)=TV_non_isotrope_separable(x5);
TV3(5)=TV_isotrope_non_separable(x5);
TV4(5)=TV_isotrope_separable(x5);

// pixel en 1 et 2
x6=x;
x6(1,1,2)=1.0;
x6(2,1,2)=1.0;
TV1(6)=TV_non_isotrope_non_separable(x6);
TV2(6)=TV_non_isotrope_separable(x6);
TV3(6)=TV_isotrope_non_separable(x6);
TV4(6)=TV_isotrope_separable(x6);

// pixel en 1 et 3
x7=x;
x7(1,1,2)=1.0;
x7(1,2,2)=1.0;
TV1(7)=TV_non_isotrope_non_separable(x7);
TV2(7)=TV_non_isotrope_separable(x7);
TV3(7)=TV_isotrope_non_separable(x7);
TV4(7)=TV_isotrope_separable(x7);


// n=3;
// x=array(double,n,n,n);
// x(1:2,1:2,3)=1.0;
// x(2,2,3)=0.0;

// ncase=11;
// TV1=array(double,ncase);
// TV2=array(double,ncase);
// TV3=array(double,ncase);
// TV4=array(double,ncase);

// // pixel en 1 => le cas que l'on veut !
// x1=x;
// x1(1,1,2)=1.0;
// TV1(1)=TV_non_isotrope_non_separable(x1);
// TV2(1)=TV_non_isotrope_separable(x1);
// TV3(1)=TV_isotrope_non_separable(x1);
// TV4(1)=TV_isotrope_separable(x1);

// // pixel en 2
// x2=x;
// x2(2,1,2)=1.0;
// TV1(2)=TV_non_isotrope_non_separable(x2);
// TV2(2)=TV_non_isotrope_separable(x2);
// TV3(2)=TV_isotrope_non_separable(x2);
// TV4(2)=TV_isotrope_separable(x2);

// // pixel en 3
// x3=x;
// x3(1,2,2)=1.0;
// TV1(3)=TV_non_isotrope_non_separable(x3);
// TV2(3)=TV_non_isotrope_separable(x3);
// TV3(3)=TV_isotrope_non_separable(x3);
// TV4(3)=TV_isotrope_separable(x3);

// // pixel en 4
// x4=x;
// x4(2,2,2)=1.0;
// TV1(4)=TV_non_isotrope_non_separable(x4);
// TV2(4)=TV_non_isotrope_separable(x4);
// TV3(4)=TV_isotrope_non_separable(x4);
// TV4(4)=TV_isotrope_separable(x4);

// // 2 pixels en 1 et 2
// x5=x;
// x5(1:2,1,2)=1.0;
// TV1(5)=TV_non_isotrope_non_separable(x5);
// TV2(5)=TV_non_isotrope_separable(x5);
// TV3(5)=TV_isotrope_non_separable(x5);
// TV4(5)=TV_isotrope_separable(x5);

// // 2 pixels en 1 et 3
// x6=x;
// x6(1,1:2,2)=1.0;
// TV1(6)=TV_non_isotrope_non_separable(x6);
// TV2(6)=TV_non_isotrope_separable(x6);
// TV3(6)=TV_isotrope_non_separable(x6);
// TV4(6)=TV_isotrope_separable(x6);

// // // 2 pixels en 3 et 4
// x7=x;
// x7(1:2,2,2)=1.0;
// TV1(7)=TV_non_isotrope_non_separable(x7);
// TV2(7)=TV_non_isotrope_separable(x7);
// TV3(7)=TV_isotrope_non_separable(x7);
// TV4(7)=TV_isotrope_separable(x7);

// // // 2 pixels en 2 et 4
// x8=x;
// x8(2,1:2,2)=1.0;
// TV1(8)=TV_non_isotrope_non_separable(x8);
// TV2(8)=TV_non_isotrope_separable(x8);
// TV3(8)=TV_isotrope_non_separable(x8);
// TV4(8)=TV_isotrope_separable(x8);

// // // 2 pixels en 1 et 4
// x9=x;
// x9(2,2,2)=1.0;
// x9(1,1,2)=1.0;
// TV1(9)=TV_non_isotrope_non_separable(x9);
// TV2(9)=TV_non_isotrope_separable(x9);
// TV3(9)=TV_isotrope_non_separable(x9);
// TV4(9)=TV_isotrope_separable(x9);

// // 3 pixels en 1, 2 et 3
// x10=x;
// x10(1:2,1:2,2)=1.0;
// x10(2,2,2)=0.0;
// TV1(10)=TV_non_isotrope_non_separable(x10);
// TV2(10)=TV_non_isotrope_separable(x10);
// TV3(10)=TV_isotrope_non_separable(x10);
// TV4(10)=TV_isotrope_separable(x10);

// // 4 pixels
// x11=x;
// x11(1:2,1:2,2)=1.0;
// TV1(11)=TV_non_isotrope_non_separable(x11);
// TV2(11)=TV_non_isotrope_separable(x11);
// TV3(11)=TV_isotrope_non_separable(x11);
// TV4(11)=TV_isotrope_separable(x11);


window, 0; fma; plg, TV1;
plg, TV2, color="blue";
plg, TV3, color="red";
plg, TV4, color="red", type="dash";

