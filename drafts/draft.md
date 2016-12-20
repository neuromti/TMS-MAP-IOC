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

In each session, TMS was applied with a MagPro X100 with MagOption and a
MCF-B70 figure-of-eight coil of 97mm diameter (Magventure, Denmark).
Stimulation was performed based on a template headmodel using the
([]{#OLE_LINK1 .anchor}TMS Navigator, Localite, Germany). We recorded
MEPs at the extensor digitorum of the right hand with surface
electromyography. The muscle was located by palpation during extension
of the wrist and anatomical landmarks. The skin was cleaned using 95%
ethanol and abrasive gel. We used self-adhesive electrodes (Neuroline,
Ambu, Germany). During the whole experiment, the subject was seated in a
comfortable reclining chair and was told to relax his or her muscles.

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

*Resting Motor Threshold *

At this spot the resting motor threshold (RMT) was determined with the
relative frequency method \[3\], i.e. RMT was defined as the stimulator
intensity at which 5 out of 10 stimuli would elicit an MEP with an
amplitude larger than 50 µVpp.

*Mapping Grid*

For mapping, a grid was created via Localite with its center 1cm
anterior of the previously detected hotspot. In steps of 0.5 cm, we set
7 x 15 grid points. This resulted in a 3 cm wide grid spanning 4.5 cm
anterior and 2.5 cm posterior to the hotspot. At each of the 105 grid
points, 3 stimuli were applied. This resulted in an average of 15
Stimuli/cm^2^. The mapping was carried out for each stimulation
condition with 110% of RMT.

*Input-Output-Curves*

After mapping, we selected the point in the grid with the highest MEP
averaged over 3 stimuli (M1 spot) and identified the most anterior point
eliciting any MEP (NPMA spot). To assess differences in the respective
neuronal structures, we measured the stimulus-response-curve at each
spot. We stimulated with 7 intensities ranging from 90% - 150% of the
respective RMT in steps of 10%. We delivered 10 stimuli per intensity.
The order of stimulation intensities was randomized within each subject.

*Signal Processing *

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
style="font-variant:small-caps;">Target</span> (M1 vs. NPMA spot) on
latency and amplitude, accounting for interactions between the three
factors and <span style="font-variant:small-caps;">Subject</span> as a
random factor. This test was performed for each stimulus-intensity. For
the time-course, we performed this analysis additionally for every
time-point. Regarding the mapping, we calculated the influence of the
categorical factors coil <span
style="font-variant:small-caps;">Orientation</span> and stimulus <span
style="font-variant:small-caps;">Waveform</span> on resting motor
threshold, latency and amplitude, accounting for interactions between
the two factors and <span
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
= -36.9, Y = -18.6. This is not different from the position of M1 as
established in literature \[4\] as evidenced by finding no significant
differences (t(12) = \[0.2, 0.39\], p = \[0.85, 0.70\]). This suggests
that, as designed, the average grid origin was 1 cm anterior to M1.

*Resting Motor Threshold*

Inspection of the average motor threshold in %MSO for biphasic
stimulation at 90° (M = 38.3, SD = 7.7) and 45° (M = 36.7, SD = 10.7),
as well as for monophasic at 90° (M = 65.9, SD = 12.1) and 45° (M =
61.7, SD = 13.0) exhibits the decreased resting motor threshold for
biphasic stimulation. Indeed, resting motor threshold was not
significantly different for orientation (F(1, 28) = 1.9, p=0.18) , but
only for waveform (F(1, 28) = 116.4, p&gt;0.001). Additionally, we found
no significant interactions (F(1, 28) = 0.1, p=0.73).

*Mapping Grid*

We found significant clusters highlighting the influence of <span
style="font-variant:small-caps;">waveform</span> and <span
style="font-variant:small-caps;">orientation</span> on MEP parameters
(figure 1). We found that biphasic instead of monophasic stimulation
increases amplitude (p = 0.047, figure 1A) when performed over a primary
motor cluster (centered on X = 33.8, Y = 23.4), while stimulation with
90° decreases amplitude (p = 0.024, figure 1D) when performed over an
anterior cluster (centered on X = 41.4, Y = 23.9). We also found that
latency was decreased during stimulation at 90° in contrast to 45° (p =
0.036, figure 1B) when performed over an anterior cluster (centered on X
= 41.1, Y = 7.7). Last, we also found evidence for an interaction
between <span style="font-variant:small-caps;">waveform</span> and <span
style="font-variant:small-caps;">orientation</span> suggesting that
biphasic at 45° and monophasic at 90° reduce latency in contrast to
biphasic at 90° and monophasic at 45° (p = 0.001, figure 1C) when
performed over primary motor areas (centered on X = -28.9, Y = -21.9).
