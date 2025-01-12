---
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
cd  /Users/driscoll/Documents/GitHub/fnc/matlab
FNC_init
```
[**Demo %s**](#demo-systems-coupledpendula)

Let's implement the coupled pendulums from {numref}`Example {number} <example-systems-coupledpendula>`. The pendulums will be pulled in opposite directions and then released together from rest.

```{literalinclude} f63_pendulums.m
:language: matlab
```

```{code-cell}
u0 = [1.25; -0.5; 0; 0];
a = 0; b = 50;
```

First we check the behavior of the system when the pendulums are uncoupled, i.e., when $k=0$.
```{tip}
:class: dropdown
Here `OutputVariables` is used to restrict output to just $u_1$ and $u_2$.
```

```{code-cell}
params =[0.01, 0.5, 0];    % gamma, L, k
ivp = ode(ODEFcn=@f63_pendulums, InitialValue=u0, Parameters=params);
theta = solutionFcn(ivp, a, b, OutputVariables = 1:2);
t = linspace(a, b, 1001);
clf, plot(t, theta(t))
xlabel("t");  ylabel("angle")
title("Uncoupled pendulums")
legend("\theta_1", "\theta_2");
```

You can see that the pendulums swing independently. Because the model is nonlinear and the initial angles are not small, they have slightly different periods of oscillation, and they go in and out of phase.

With coupling activated, a different behavior is seen.

```{code-cell}
params(3) = 1;
ivp = ode(ODEFcn=@f63_pendulums, InitialValue=u0, Parameters=params);
theta = solutionFcn(ivp, a, b, OutputVariables = 1:2);
clf, plot(t, theta(t))
xlabel("t");  ylabel("angle")
title("Coupled pendulums")
legend("\theta_1", "\theta_2");
```

The coupling makes the pendulums swap energy back and forth.
