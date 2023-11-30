within ;
model EX2_2Switch
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor heatCapacitor(C=3000, T(
        displayUnit="degC", start=393.15))
    annotation (Placement(transformation(extent={{-10,100},{10,120}})));
  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow fixedHeatFlow(
    Q_flow=2414.5,
    T_ref=293.15,
    alpha=0)
    annotation (Placement(transformation(extent={{-92,68},{-72,88}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G=
        100) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-38,78})));
  Modelica.Thermal.FluidHeatFlow.Components.OpenTank openTank(
    medium=Modelica.Thermal.FluidHeatFlow.Media.Water_10degC(),
    T0(displayUnit="degC") = 283.15,
    T0fixed=true,
    ATank=0.01,
    hTank=0.8,
    pAmbient=100000,
    useHeatPort=true,
    level(start=0.8, fixed=true))
    annotation (Placement(transformation(extent={{-242,26},{-222,46}})));
  Modelica.Thermal.FluidHeatFlow.Components.Pipe pipe(
    medium=Modelica.Thermal.FluidHeatFlow.Media.Water_10degC(),
    m=0.502,
    T0(displayUnit="K") = 283.15,
    T0fixed=false,
    dpLaminar=0,
    dpNominal=0,
    frictionLoss=0,
    useHeatPort=true,
    h_g=0) annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=180)));
  Modelica.Thermal.HeatTransfer.Components.Convection convection
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=270,
        origin={0,50})));
  Modelica.Blocks.Sources.Constant const(k=300) annotation (Placement(
        transformation(
        extent={{10,-10},{-10,10}},
        rotation=0,
        origin={36,50})));
  Modelica.Thermal.FluidHeatFlow.Sources.VolumeFlow volumeFlow(
    medium=Modelica.Thermal.FluidHeatFlow.Media.Water_10degC(),
    T0(displayUnit="K") = 283.15,
    V_flow(start=0),
    useVolumeFlowInput=true,
    constantVolumeFlow=0) annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=0,
        origin={-116,0})));
  Modelica.Thermal.HeatTransfer.Celsius.TemperatureSensor
                                         temperatureSensor annotation (
      Placement(transformation(
        origin={26,78},
        extent={{-10,-10},{10,10}},
        rotation=0)));
  Modelica.Blocks.Logical.Hysteresis hysteresisT(
    uLow=42,
    uHigh=59,
    pre_y_start=false)
    annotation (Placement(transformation(extent={{56,68},{76,88}})));
  Modelica.Blocks.Logical.Switch switch1 annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=180,
        origin={-24,-44})));
  Modelica.Blocks.Sources.Constant const1(k=5e-5) annotation (Placement(
        transformation(
        extent={{10,-10},{-10,10}},
        rotation=0,
        origin={42,-24})));
  Modelica.Blocks.Sources.Constant const2(k=0) annotation (Placement(
        transformation(
        extent={{10,-10},{-10,10}},
        rotation=0,
        origin={42,-62})));
  Modelica.Thermal.FluidHeatFlow.Components.OpenTank openTank2(
    medium=Modelica.Thermal.FluidHeatFlow.Media.Water_10degC(),
    T0(displayUnit="K") = 283.15,
    T0fixed=true,
    ATank=0.01,
    hTank=0.8,
    pAmbient=100000,
    useHeatPort=true,
    level(start=1e-12, fixed=true))
    annotation (Placement(transformation(extent={{120,26},{140,46}})));
  Modelica.Blocks.Logical.Hysteresis hysteresisT1(uLow=273.15 + 5, uHigh=273.15
         + 10)
    annotation (Placement(transformation(extent={{-206,38},{-186,58}})));
  Modelica.Blocks.Logical.Switch switch2 annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=180,
        origin={-188,76})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-266,76})));
  Modelica.Blocks.Sources.Constant const3(k=-4000) annotation (Placement(
        transformation(
        extent={{10,-10},{-10,10}},
        rotation=0,
        origin={-126,96})));
  Modelica.Blocks.Sources.Constant const4(k=0) annotation (Placement(
        transformation(
        extent={{10,-10},{-10,10}},
        rotation=0,
        origin={-116,38})));
  Modelica.Blocks.Logical.Hysteresis hysteresisT2(uLow=5.5 + 273.15, uHigh=9 +
        273.15)
    annotation (Placement(transformation(extent={{202,36},{222,56}})));
  Modelica.Blocks.Logical.Switch switch3 annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=180,
        origin={220,74})));
  Modelica.Blocks.Sources.Constant const5(k=-4000) annotation (Placement(
        transformation(
        extent={{10,-10},{-10,10}},
        rotation=0,
        origin={282,94})));
  Modelica.Blocks.Sources.Constant const6(k=0) annotation (Placement(
        transformation(
        extent={{10,-10},{-10,10}},
        rotation=0,
        origin={292,36})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow prescribedHeatFlow1
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={174,74})));
  Modelica.Thermal.FluidHeatFlow.Components.Valve F2(
    medium=Modelica.Thermal.FluidHeatFlow.Media.Water_10degC(),
    T0=283.15,
    T0fixed=true,
    y1=1,
    Kv1=1,
    dp0=0) annotation (Placement(transformation(extent={{-86,10},{-66,-10}})));
  Modelica.Thermal.FluidHeatFlow.Components.Valve B1(
    medium=Modelica.Thermal.FluidHeatFlow.Media.Water_10degC(),
    T0=283.15,
    T0fixed=true,
    y1=1,
    Kv1=1,
    dp0=0)
    annotation (Placement(transformation(extent={{-120,-54},{-100,-74}})));
  Modelica.Thermal.FluidHeatFlow.Components.Valve B2(
    medium=Modelica.Thermal.FluidHeatFlow.Media.Water_10degC(),
    T0=283.15,
    T0fixed=true,
    y1=1,
    Kv1=1,
    dp0=0) "B2"
    annotation (Placement(transformation(extent={{-186,28},{-166,8}})));
  Modelica.Thermal.FluidHeatFlow.Components.Valve F1(
    medium=Modelica.Thermal.FluidHeatFlow.Media.Water_10degC(),
    T0=283.15,
    T0fixed=true,
    y1=1,
    Kv1=1,
    dp0=0) annotation (Placement(transformation(extent={{-226,10},{-206,-10}})));
  Modelica.Blocks.Logical.Hysteresis hysteresisSwitchTank(
    uLow=1e-12,
    uHigh=0.8,
    pre_y_start=true)
    annotation (Placement(transformation(extent={{-346,-70},{-326,-50}})));
  Modelica.Blocks.Logical.Not not1 annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-294,-116})));
  Modelica.Blocks.Math.BooleanToReal booleanToReal
    annotation (Placement(transformation(extent={{-274,-70},{-254,-50}})));
  Modelica.Blocks.Math.BooleanToReal booleanToReal1
    annotation (Placement(transformation(extent={{-274,-126},{-254,-106}})));
equation
  connect(heatCapacitor.port, thermalConductor.port_b)
    annotation (Line(points={{0,100},{0,78},{-28,78}}, color={191,0,0}));
  connect(fixedHeatFlow.port, thermalConductor.port_a)
    annotation (Line(points={{-72,78},{-48,78}}, color={191,0,0}));
  connect(convection.fluid, pipe.heatPort) annotation (Line(points={{
          -1.77636e-15,40},{-1.77636e-15,24},{1.72085e-15,24},{1.72085e-15,10}},
        color={191,0,0}));
  connect(const.y, convection.Gc)
    annotation (Line(points={{25,50},{10,50}}, color={0,0,127}));
  connect(convection.solid, thermalConductor.port_b) annotation (Line(points={{
          1.77636e-15,60},{0,60},{0,78},{-28,78}}, color={191,0,0}));
  connect(temperatureSensor.port, thermalConductor.port_b)
    annotation (Line(points={{16,78},{-28,78}}, color={191,0,0}));
  connect(hysteresisT.u, temperatureSensor.T)
    annotation (Line(points={{54,78},{36,78}}, color={0,0,127}));
  connect(hysteresisT.y, switch1.u2) annotation (Line(points={{77,78},{82,78},{
          82,-44},{-12,-44}}, color={255,0,255}));
  connect(const1.y, switch1.u1) annotation (Line(points={{31,-24},{-6,-24},{-6,
          -36},{-12,-36}}, color={0,0,127}));
  connect(const2.y, switch1.u3) annotation (Line(points={{31,-62},{-4,-62},{-4,
          -52},{-12,-52}}, color={0,0,127}));
  connect(switch1.y, volumeFlow.volumeFlow) annotation (Line(points={{-35,-44},
          {-116,-44},{-116,-10}}, color={0,0,127}));
  connect(pipe.flowPort_a, openTank2.flowPort) annotation (Line(points={{10,
          -1.77636e-15},{130,-1.77636e-15},{130,26}}, color={255,0,0}));
  connect(hysteresisT1.y, switch2.u2) annotation (Line(points={{-185,48},{-134,
          48},{-134,76},{-176,76}}, color={255,0,255}));
  connect(prescribedHeatFlow.Q_flow, switch2.y)
    annotation (Line(points={{-256,76},{-199,76}}, color={0,0,127}));
  connect(const4.y, switch2.u3) annotation (Line(points={{-127,38},{-130,38},{
          -130,52},{-170,52},{-170,68},{-176,68}}, color={0,0,127}));
  connect(const3.y, switch2.u1) annotation (Line(points={{-137,96},{-170,96},{
          -170,84},{-176,84}}, color={0,0,127}));
  connect(prescribedHeatFlow.port, openTank.heatPort) annotation (Line(points={
          {-276,76},{-282,76},{-282,26},{-242,26}}, color={191,0,0}));
  connect(hysteresisT2.y, switch3.u2) annotation (Line(points={{223,46},{274,46},
          {274,74},{232,74}}, color={255,0,255}));
  connect(const6.y, switch3.u3) annotation (Line(points={{281,36},{278,36},{278,
          50},{238,50},{238,66},{232,66}}, color={0,0,127}));
  connect(const5.y, switch3.u1) annotation (Line(points={{271,94},{238,94},{238,
          82},{232,82}}, color={0,0,127}));
  connect(hysteresisT2.u, openTank2.TTank) annotation (Line(points={{200,46},{
          146,46},{146,30},{141,30}}, color={0,0,127}));
  connect(switch3.y, prescribedHeatFlow1.Q_flow)
    annotation (Line(points={{209,74},{184,74}}, color={0,0,127}));
  connect(prescribedHeatFlow1.port, openTank2.heatPort) annotation (Line(points
        ={{164,74},{106,74},{106,26},{120,26}}, color={191,0,0}));
  connect(hysteresisT1.u, openTank.TTank) annotation (Line(points={{-208,48},{
          -216,48},{-216,30},{-221,30}}, color={0,0,127}));
  connect(volumeFlow.flowPort_b, F2.flowPort_a)
    annotation (Line(points={{-106,0},{-86,0}}, color={255,0,0}));
  connect(F2.flowPort_b, pipe.flowPort_b) annotation (Line(points={{-66,0},{-42,
          0},{-42,7.21645e-16},{-10,7.21645e-16}}, color={255,0,0}));
  connect(F2.flowPort_b, B1.flowPort_b) annotation (Line(points={{-66,0},{-52,0},
          {-52,-64},{-100,-64}}, color={255,0,0}));
  connect(B2.flowPort_b, F2.flowPort_a) annotation (Line(points={{-166,18},{-96,
          18},{-96,0},{-86,0}}, color={255,0,0}));
  connect(B2.flowPort_a, openTank.flowPort)
    annotation (Line(points={{-186,18},{-232,18},{-232,26}}, color={255,0,0}));
  connect(F1.flowPort_a, openTank.flowPort)
    annotation (Line(points={{-226,0},{-232,0},{-232,26}}, color={255,0,0}));
  connect(F1.flowPort_b, B1.flowPort_a) annotation (Line(points={{-206,0},{-156,
          0},{-156,-64},{-120,-64}}, color={255,0,0}));
  connect(F1.flowPort_b, volumeFlow.flowPort_a)
    annotation (Line(points={{-206,0},{-126,0}}, color={255,0,0}));
  connect(openTank.level, hysteresisSwitchTank.u) annotation (Line(points={{
          -221,36},{-220,36},{-220,60},{-390,60},{-390,-60},{-348,-60}}, color=
          {0,0,127}));
  connect(booleanToReal.u, hysteresisSwitchTank.y)
    annotation (Line(points={{-276,-60},{-325,-60}}, color={255,0,255}));
  connect(not1.u, hysteresisSwitchTank.y) annotation (Line(points={{-306,-116},
          {-314,-116},{-314,-60},{-325,-60}}, color={255,0,255}));
  connect(not1.y, booleanToReal1.u)
    annotation (Line(points={{-283,-116},{-276,-116}}, color={255,0,255}));
  connect(booleanToReal1.y, B1.y) annotation (Line(points={{-253,-116},{-110,
          -116},{-110,-74}}, color={0,0,127}));
  connect(B2.y, booleanToReal1.y) annotation (Line(points={{-176,8},{-176,-116},
          {-253,-116}}, color={0,0,127}));
  connect(booleanToReal.y, F1.y) annotation (Line(points={{-253,-60},{-216,-60},
          {-216,-10}}, color={0,0,127}));
  connect(F2.y, booleanToReal.y) annotation (Line(points={{-76,-10},{-76,-86},{
          -216,-86},{-216,-60},{-253,-60}}, color={0,0,127}));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-280,-100},{280,
            180}})),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-280,-100},{
            280,180}}), graphics={Text(
          extent={{-18,-10},{28,-24}},
          textColor={28,108,200},
          textString="migliora massa
")}),
    uses(Modelica(version="4.0.0")));
end EX2_2Switch;
