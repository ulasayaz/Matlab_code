function draw2D(x);

%%% plots the 

s=length(x);
if s == 2
    X=zeros(2,1);
    Y=zeros(2,1);

    X(1)=x(1);
    Y(1)=x(2);

    plot(X,Y) % can also use 'line' function
elseif s == 1
    X=zeros(2,1);
    Y=zeros(2,1);

    X(1)=real(x);
    Y(1)=imag(x);

    plot(X,Y) % can also use 'line' function

else
    'dim more than 2'
end

