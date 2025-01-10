using HDF5,Plots

default(dpi=400,linewidth=1.5,label="")
W1 = vec(h5read("stabreg_data_bd.h5","W1"));
W2 = vec(h5read("stabreg_data_bd.h5","W2"));
W3 = vec(h5read("stabreg_data_bd.h5","W3"));
W4 = vec(h5read("stabreg_data_bd.h5","W4"));

plt = plot(aspect_ratio=1,size=(760,380),frame=:zerolines,layout=(1,2),
	xaxis=((-4,12),-8:4:12),yaxis=((-8,8),-8:4:8) )
plot!(Shape([20; real(W1); 20; 20; -20; -20],[20; imag(W1); 20; -20; -20; 20]),fillcolor=RGBA(.7,.7,1,.25),l=0,subplot=1)
plot!(Shape([20; real(W2); 20; 20; -20; -20],[20; imag(W2); 20; -20; -20; 20]),fillcolor=RGBA(.7,.7,1,.25),l=0,subplot=1)
plot!(Shape([20; real(W3); 20; 20; -20; -20],[20; imag(W3); 20; -20; -20; 20]),fillcolor=RGBA(.7,.7,1,.25),l=0,subplot=1)
plot!(Shape([20; real(W4); 20; 20; -20; -20],[20; imag(W4); 20; -20; -20; 20]),fillcolor=RGBA(.7,.7,1,.25),l=0,subplot=1)
plot!(real(W1),imag(W1),linecolor=:darkblue,subplot=1);
plot!(real(W2),imag(W2),linecolor=:darkblue,subplot=1);
plot!(real(W3),imag(W3),linecolor=:darkblue,subplot=1);
plot!(real(W4),imag(W4),linecolor=:darkblue,subplot=1);

annotate!(2.2,.73,text("1",:bold,9),subplot=1);
annotate!(4.2,1.4,text("2",:bold,9),subplot=1);
annotate!(6.5,2.4,text("3",:bold,9),subplot=1);
annotate!(10.1,3.95,text("4",:bold,9),subplot=1)

W1 = vec(h5read("stabreg_data_rk.h5","W1"));
W2 = vec(h5read("stabreg_data_rk.h5","W2"));
W3 = vec(h5read("stabreg_data_rk.h5","W3"));
W4 = vec(h5read("stabreg_data_rk.h5","W4"));

plot!(subplot=2,xaxis=((-4,2),-4:2),yaxis=((-3,3),-3:3))
plot!(Shape(real(W1),imag(W1)),fillcolor=RGBA(.7,.7,1,.25),linecolor=:darkblue,subplot=2);
plot!(Shape(real(W2),imag(W2)),fillcolor=RGBA(.7,.7,1,.25),linecolor=:darkblue,subplot=2);
plot!(Shape(real(W3),imag(W3)),fillcolor=RGBA(.7,.7,1,.25),linecolor=:darkblue,subplot=2);
plot!(Shape(real(W4),imag(W4)),fillcolor=RGBA(.7,.7,1,.25),linecolor=:darkblue,subplot=2);

annotate!(-.51,.65,text("1",:bold,9),subplot=2);
annotate!(-.45,1.28,text("2",:bold,9),subplot=2);
annotate!(-.3,2,text("3",:bold,9),subplot=2);
annotate!(-.16,2.75,text("4",:bold,9),subplot=2)

savefig(plt,"stabreg_bd_rk.pdf")
savefig(plt,"stabreg_bd_rk.svg")