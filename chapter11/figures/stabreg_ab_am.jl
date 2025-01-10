using HDF5,Plots

default(dpi=400,linewidth=1.5,label="")
W1 = vec(h5read("stabreg_data_ab.h5","W1"));
W2 = vec(h5read("stabreg_data_ab.h5","W2"));
W3 = vec(h5read("stabreg_data_ab.h5","W3"));
W4 = vec(h5read("stabreg_data_ab.h5","W4"));

plt = plot(aspect_ratio=1,size=(760,380),frame=:zerolines,layout=(1,2),
	xaxis=((-2.25,0.25),-3:0.5:1),yaxis=((-1.25,1.25),-2:.5:2) )
plot!(Shape(real(W1),imag(W1)),fillcolor=RGBA(.7,.7,1,.25),linecolor=:darkblue,subplot=1);
plot!(Shape(real(W2),imag(W2)),fillcolor=RGBA(.7,.7,1,.25),linecolor=:darkblue,subplot=1);
plot!(Shape(real(W3),imag(W3)),fillcolor=RGBA(.7,.7,1,.25),linecolor=:darkblue,subplot=1);
plot!(Shape(real(W4),imag(W4)),fillcolor=RGBA(.7,.7,1,.25),linecolor=:darkblue,subplot=1);

annotate!(-1.5,.77,text("1",:bold,9),subplot=1);
annotate!(-.7,.55,text("2",:bold,9),subplot=1);
annotate!(-.36,.36,text("3",:bold,9),subplot=1);
annotate!(-.16,.22,text("4",:bold,9),subplot=1);


W3 = vec(h5read("stabreg_data_am.h5","W3"));
W4 = vec(h5read("stabreg_data_am.h5","W4"));
W5 = vec(h5read("stabreg_data_am.h5","W5"));

plot!(subplot=2,xaxis=((-7,1),-6:2:2),yaxis=((-4,4),-4:2:4))
plot!(Shape([0, -10, -10, 0],[5, 5, -5, -5]),fillcolor=RGBA(.7,.7,1,.25),linecolor=:darkblue,subplot=2)
plot!(Shape(real(W3),imag(W3)),fillcolor=RGBA(.7,.7,1,.25),linecolor=:darkblue,subplot=2);
plot!(Shape(real(W4),imag(W4)),fillcolor=RGBA(.7,.7,1,.25),linecolor=:darkblue,subplot=2);
plot!(Shape(real(W5),imag(W5)),fillcolor=RGBA(.7,.7,1,.25),linecolor=:darkblue,subplot=2);

annotate!(-.3,3.7,text("2",:bold,9),subplot=2);
annotate!(-4.6,2.3,text("3",:bold,9),subplot=2);
annotate!(-2.5,1,text("4",:bold,9),subplot=2);
annotate!(-1.4,.75,text("5",:bold,9),subplot=2);

savefig(plt,"stabreg_ab_am.pdf")
savefig(plt,"stabreg_ab_am.svg")