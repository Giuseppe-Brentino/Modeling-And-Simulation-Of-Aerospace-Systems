within ;
model EX_2_1
  Modelica.Electrical.Analog.Basic.Resistor resistor(R=0.1) annotation (
      Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-48,38})));
  Modelica.Electrical.Analog.Basic.Inductor inductor(L=0.01) annotation (
      Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-20,38})));
  Modelica.Electrical.Analog.Basic.RotationalEMF emf(k=0.3) annotation (
      Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-2,22})));
  Modelica.Mechanics.Rotational.Components.Inertia inertia(J=0.001)
    annotation (Placement(transformation(extent={{22,12},{42,32}})));
  Modelica.Electrical.Analog.Basic.Ground ground annotation (Placement(
        transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-2,-22})));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G=
        100) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={54,-8})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor(C=3000)
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={88,-26})));
  Modelica.Mechanics.Rotational.Components.Inertia Propeller(
    J=1/12*(0.8^3)*2700*Modelica.Constants.pi*(0.01/2)^2,
    phi(start=0),
    w(start=0),
    a(start=0))
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={162,22})));
  Modelica.Blocks.Continuous.PI PI(
    k=0.07,
    T=0.05,
    initType=Modelica.Blocks.Types.Init.InitialState,
    x_start=0)                                   annotation (Placement(
        transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-146,10})));
  Modelica.Blocks.Sources.Step step(height=210, startTime=5) annotation (
      Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-242,10})));
  Modelica.Blocks.Math.Feedback feedback annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-202,10})));
  Modelica.Mechanics.Rotational.Sensors.SpeedSensor speedSensor annotation (
      Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-24,-58})));
  Modelica.Electrical.Analog.Sources.SignalVoltage signalVoltage annotation (
      Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={-96,10})));
  Modelica.Mechanics.Rotational.Sources.QuadraticSpeedDependentTorque
    quadraticSpeedDependentTorque(
    tau_nominal=-100,
    TorqueDirection=false,                         w_nominal=210) annotation (
      Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={220,22})));
  Modelica.Mechanics.Rotational.Components.LossyGear lossyGear(
    useSupport=false,
    ratio=2,
    lossTable=[0,0.99,0.99,0,0; 50,0.98,0.98,0.5,0.5; 100,0.97,0.97,1,1; 210,
        0.96,0.96,1.5,1.5],
    useHeatPort=true)
    annotation (Placement(transformation(extent={{54,12},{74,32}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
    annotation (Placement(transformation(extent={{108,-6},{128,14}})));
equation
  connect(emf.flange, inertia.flange_a)
    annotation (Line(points={{8,22},{10,22},{10,22},{22,22}}, color={0,0,0}));
  connect(ground.p, emf.n)
    annotation (Line(points={{-2,-12},{-2,12}}, color={0,0,255}));
  connect(inductor.p, resistor.n)
    annotation (Line(points={{-30,38},{-38,38}}, color={0,0,255}));
  connect(inductor.n, emf.p)
    annotation (Line(points={{-10,38},{-2,38},{-2,32}}, color={0,0,255}));
  connect(speedSensor.w, feedback.u2) annotation (Line(points={{-35,-58},{-202,
          -58},{-202,2}},  color={0,0,127}));
  connect(step.y, feedback.u1)
    annotation (Line(points={{-231,10},{-210,10}}, color={0,0,127}));
  connect(feedback.y, PI.u)
    annotation (Line(points={{-193,10},{-158,10}}, color={0,0,127}));
  connect(speedSensor.flange, Propeller.flange_b) annotation (Line(points={{-14,-58},
          {178,-58},{178,22},{172,22}},      color={0,0,0}));
  connect(signalVoltage.n, ground.p) annotation (Line(points={{-96,0},{-96,-6},
          {-2,-6},{-2,-12}},               color={0,0,255}));
  connect(PI.y, signalVoltage.v) annotation (Line(points={{-135,10},{-108,10}},
                                     color={0,0,127}));
  connect(resistor.p, signalVoltage.p) annotation (Line(points={{-58,38},{-96,
          38},{-96,20}},          color={0,0,255}));
  connect(thermalConductor.port_b, heatCapacitor.port)
    annotation (Line(points={{54,-18},{54,-26},{78,-26}}, color={191,0,0}));
  connect(quadraticSpeedDependentTorque.flange, Propeller.flange_b)
    annotation (Line(points={{210,22},{172,22}}, color={0,0,0}));
  connect(lossyGear.flange_b, Propeller.flange_a)
    annotation (Line(points={{74,22},{152,22}},color={0,0,0}));
  connect(lossyGear.flange_a, inertia.flange_b)
    annotation (Line(points={{54,22},{42,22}}, color={0,0,0}));
  connect(lossyGear.heatPort, thermalConductor.port_a)
    annotation (Line(points={{54,12},{54,2}}, color={191,0,0}));
  connect(temperatureSensor.port, thermalConductor.port_a)
    annotation (Line(points={{108,4},{70,4},{70,2},{54,2}}, color={191,0,0}));
  annotation (uses(Modelica(version="4.0.0")), Diagram(graphics={Text(extent={{
              94,58},{94,58}}, textColor={28,108,200})}));
end EX_2_1;
