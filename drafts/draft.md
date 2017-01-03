**Materials and Methods**

[]{#OLE_LINK2 .anchor}*Study Design*

13 healthy subjects participated in this experiment (age M = 25.7, SD =
7.1 years; 7 male). All subjects were right-handed and scored at least
+75% on the Edinburgh Handedness Inventory \[1\]. None of them had any
history of neurological or psychiatric illnesses or had any
contraindications to TMS \[2\]. Subjects gave written informed consent
before participation, and the study was approved by the local ethics
committee.

We performed hotspot detection and a cortical mapping for each of four
stimulation conditions: biphasic latero-medial (BI-90), biphasic
posterolateral-anteromedial (BI-45), monophasic latero-medial (MO-90),
and monophasic posterolateral-anteromedial (MO-45). Due to the duration
of each mapping, recordings were performed in two sessions. Per session,
two conditions were assessed. The order of conditions was randomized
within each subject.

*Transcranial Magnetic Stimulation*

In each session, TMS was applied with a MagPro X100 with MagOption using
Power Mode and a MCF-B70 figure-of-eight coil of 97mm diameter
(Magventure, Denmark). Consider that Power Mode increases the maximum
stimulator output by 40%. Stimulation was performed based on a template
headmodel using the ([]{#OLE_LINK1 .anchor}TMS Navigator, Localite,
Germany). We recorded MEPs at the extensor digitorum of the right hand
with surface electromyography. The muscle was located by palpation
during extension of the wrist and anatomical landmarks. The skin was
cleaned using 95% ethanol and abrasive gel. We used self-adhesive
electrodes (Neuroline, Ambu, Germany). During the whole experiment, the
subject was seated in a comfortable reclining chair and was told to
relax his or her muscles.

*Hotspot Detection*

The individual hotspot for EDC representation in the right primary motor
cortex was determined for each of the four conditions prior to every
mapping. Based on the individual´s head anatomy, we started stimulation
at the approximate location of the hotspot. Initial stimulation
intensity was set to 40% of maximum stimulator output (MSO) for biphasic
stimulation and 65% of MSO for monophasic stimulation \[2\]. If no MEPs
could be elicited, intensity was increased in 5% steps. Coil position
was adapted in a random fashion, while coil orientation was fixed based
on the respective condition. Finally, the hotspot of each condition was
defined as the spot eliciting the highest MEP with shortest latency for
a given stimulation intensity.

*Resting Motor Threshold*

At this spot the resting motor threshold (RMT) was determined
[]{#OLE_LINK3 .anchor}with the relative frequency method \[3\], i.e. RMT
was defined as the stimulator intensity at which 5 out of 10 stimuli
would elicit an MEP with an amplitude larger than 50 µVpp.

*Mapping Grid*

For mapping, a grid was created via Localite with its center 1cm
anterior of the previously detected hotspot. In steps of 0.5 cm, we set
7 x 15 grid points. This resulted in a 3 cm wide grid spanning 4.5 cm
anterior and 2.5 cm posterior to the hotspot. At each of the 105 grid
points, 3 stimuli were applied. This resulted in an average of 15
Stimuli/cm^2^. The mapping was carried out for each stimulation
condition with 110% of the RMT of the primary target. We had to remove 6
measurements for technical reasons, resulting in 46 total and 11 to 13
maps for each of the four conditions.

*Input-Output-Curves*

After mapping, we selected the point in the grid with the highest MEP
averaged over 3 stimuli (primary target) and identified the most
anterior point eliciting any MEP (anterior target). To assess
differences in the respective neuronal structures, we measured the
stimulus-response-curve at each target. We stimulated with 7 intensities
ranging from 90% - 150% of the respective RMT in steps of 10%,
delivering 10 stimuli per intensity. The order of stimulation
intensities was randomized within each subject. Consider that 150% of
67% MSO is 100%, suggesting that at RMT higher than 67%, we would not be
able to sample for high intensities. To be able to adequately sample the
complete input-output curve, we therefore used the following procedure.
If the RMT at the primary target was above 67% MSO, we calculated the
stimulation intensities for primary and anterior target relative to 67%
MSO. If RMT at the primary target was below 67%, we determined the RMT
for the anterior target with the relative frequency method \[3\], and
used distinct stimulation intensities for primary and anterior target.
We had to remove 15 measurements for technical reasons, resulting in 89
total and 10 to 13 input-output curve measurements for each of the eight
conditions.

*Signal Processing*

For all MEPs, latency and amplitude were estimated automatically
offline, using costum-written MatLab functions. Additionally, all trials
were visually inspected for correctness of the estimation. Incorrect
estimates were corrected, and artifacted trials were rejected. Latency
and amplitude were averaged for each grid point and for each intensity
level. For input-output-curve assessment, we also processed the
time-course of MEPs from 5 to 60ms after the TMS pulse. The raw
time-course was detrended and baselined in reference to the period from
5 to 17 ms after the TMS pulse.

*Statistical Analysis*

Regarding the input-output-curve, we calculated the influence of the
categorical factors coil <span
style="font-variant:small-caps;">Orientation</span> (90° vs. 45°) and
stimulus <span style="font-variant:small-caps;">Waveform</span>
(biphasic vs. monophasic) as well as <span
style="font-variant:small-caps;">Target</span> (M1 vs. NPMA) on latency
and amplitude, accounting for interactions between the three factors and
<span style="font-variant:small-caps;">Subject</span> as a random
factor. This test was performed for each stimulus-intensity. For the
time-course, we performed this analysis additionally for every
time-point. We also calculated the influence of <span
style="font-variant:small-caps;">Target</span> on measured RMT between
anterior and primary target using an analysis of variance, accounting
for <span style="font-variant:small-caps;">Subject</span> as a random
factor.

Regarding the mapping, we calculated the influence of the categorical
factors coil <span style="font-variant:small-caps;">Orientation</span>
and stimulus <span style="font-variant:small-caps;">Waveform</span> on
resting motor threshold, latency and amplitude, accounting for
interactions between the two factors and <span
style="font-variant:small-caps;">Subject</span> as a random factor. This
test was performed for each grid-point. The statistical significance of
the influence of the factors on latency and amplitude was additionally
estimated by contrasting the model coefficients with a permutation
analysis using 1000 repetitions. For assessment of the topology, we
additionally performed a cluster-based permutation analysis based on the
sum of coefficients of neighboring significant grid points. The
significance threshold was set to 5% for all statistical tests.

**Results**

*Hotspot Detection*

The average position of the hotspot across subjects used for the
definition of the two-dimensional mapping grid origin was centered on X
= -36.9, Y = -18.6. We found no significant differences (t(12) = \[0.19,
0.39\], p = \[0.85, 0.70\]) to the position from M1 as established in
literature \[4\]. This suggests that, as designed, the average grid
origin was 1 cm anterior to M1.

*Mapping Grid*

The average stimulation intensity used for mapping was significantly
different for <span style="font-variant:small-caps;">waveform</span>
(F(1, 28) = 116.4, p&gt;0.001) but not for <span
style="font-variant:small-caps;">orientation</span> (F(1, 28) = 1.9,
p=0.18). Additionally, we found no significant interactions (F(1, 28) =
0.1, p=0.73). We found grid points exhibiting a significant influence of
<span style="font-variant:small-caps;">waveform</span> and <span
style="font-variant:small-caps;">orientation</span> on MEP parameters
(figure 1). Cluster-based permutation test revealed a significant
difference for biphasic instead of monophasic stimulation increases
amplitude (p = 0.047, figure 1A), mainly pronounced over primary motor
areas (centered on X = -33.8, Y = -23.4), while stimulation with 90°
decreases amplitude (p = 0.024, figure 1D) over anterior areas (centered
on X = -36.6, Y = 3.2). The cluster-based permutation tests also
revealed that latency was decreased during stimulation at 90° in
contrast to 45° (p = 0.036, figure 1B), mainly over anterior areas
(centered on X = -41.1, Y = 7.7). Last, we also found evidence for an
interaction between <span
style="font-variant:small-caps;">waveform</span> and <span
style="font-variant:small-caps;">orientation</span> suggesting that
biphasic at 45° and monophasic at 90° reduce latency in contrast to
biphasic at 90° and monophasic at 45° (p = 0.001, figure 1C) mainly over
primary motor areas (centered on X = -28.8, Y = -22.0).

*Anterior Stimulation Target*

The distribution of the anterior hotspot across subjects is clustered
around antero-medial grid points (see figure 2A). We projected the
positions also unto a template cortical surface for visual
representation (see figure 2B). Anterior stimulation was delivered in
average 21.1 mm anterior (CI95% = 18.8 - 23.5 mm) and 4.3 mm medial
(CI95% = 1.4 - 7.2 mm) to the M1 \[4\], with a center of gravity at X =
-32.4 and Y = 3.5. The anterior stimulation targets closest classically
motor-related area is the dorsal premotor area \[4\] (see figure 2C),
yet it is still 4.7 mm anterior (CI95% = 2.3 - 7.0 mm) to this region.

*Input-Output-Curves*

Based on the permutation test, we found no evidence for an interaction
between <span style="font-variant:small-caps;">Orientation</span>, <span
style="font-variant:small-caps;">waveform</span> and <span
style="font-variant:small-caps;">target</span> at any stimulation
intensity (all p &gt; 0.062), suggesting similarity of the response
curves for anterior and primary targets. We did find significant
modulation for the factors <span
style="font-variant:small-caps;">Orientation</span>, <span
style="font-variant:small-caps;">waveform</span> and <span
style="font-variant:small-caps;">Target</span> at a local 5% alpha error
level, but after applying Bonferonni-correction, effects were limited to
<span style="font-variant:small-caps;">waveform</span>. Biphasic
stimulation exhibited increased amplitude (see figure 3B) and decreased
latency (see figure 3E) at moderate stimulation intensities (around 120
- 130 % RMT).

![figure1.tif](media/image1.tiff){width="6.299212598425197in"
height="5.225859580052494in"}**Figure 1:** *It shows the topography of
significance for the two factors <span
style="font-variant:small-caps;">Orientation</span> and <span
style="font-variant:small-caps;">Waveform</span> on latency and
amplitude at the 105 grid points. The maps were interpolated, and colors
indicate significance level (red increased, blue decreased). Grey
contour lines indicate the threshold for significance at the 5% level.
Additionally, we added to each significant cluster a textbox with the
estimation of its p-value based on the results of the cluster
permutation analysis. **A** shows the influence of <span
style="font-variant:small-caps;">Waveform</span> amplitude, suggesting
increased amplitude for biphasic in contrast to monophasic stimulation.
**B** shows the influence of <span
style="font-variant:small-caps;">Orientation</span> on latency,
suggesting decreased latency for stimulation with 90° in contrast to
45°. **C** shows the influence of <span
style="font-variant:small-caps;">Waveform x Orientation</span> on
latency, suggesting decreased latency for stimulation with biphasic 90°
and monophasic 45° in contrast to biphasic 45° and monophasic 90°. **D**
shows the influence of <span
style="font-variant:small-caps;">Orientation</span> on amplitude,
suggesting decreased amplitude for stimulation with 90° in contrast to
45°.*

![figure2.tif](media/image2.tiff){width="6.299807524059492in"
height="3.567361111111111in"}

**Figure 2:** *It shows the spatial distribution of the anterior
stimulation point on the grid (**A**) as well as projected unto a
template headmodel (**B**). Colors indicate the probability density
estimate. **C** shows the distance of each subject's anterior
stimulation position to several motor-related brain regions with a
scatter-cloud and the mean presented as a bar.*

*\
*

![image3.tif](media/image3.tiff){width="6.299212598425197in"
height="5.640583989501312in"}

**Figure 3:** *It shows the influence of the factors <span
style="font-variant:small-caps;">Orientation</span>, <span
style="font-variant:small-caps;">Waveform</span> and <span
style="font-variant:small-caps;">Target</span> on latency and amplitude
measures for different stimulation intensities. The left column
**(A-C)** shows amplitude with µVpp on the y-axis, and the right column
**(D-F)** shows latency with ms on the y-axis. Data is shown with
stimulation intensity on the x-axis in percent of resting motor
threshold. Enlarged color-markers indicate significant differences
between the respective factor levels as returned by permutation analysis
(Bonferonni-corrected for 7 x 6 multiple comparison). Colored patches
indicate the bootstrapped 95% confidence intervals. *
