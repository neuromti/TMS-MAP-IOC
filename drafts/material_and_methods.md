**Materials and Methods**

_Study Design_

13 healthy subjects participated in this experiment (age M = 25.7, SD = 7.1 years; 7 male). All subjects were right-handed and scored at least +75% on the Edinburgh Handedness Inventory [1]. None of them had any history of neurological or psychiatric illnesses or had any contraindications to TMS [2]. Subjects gave written informed consent before participation, and the study was approved by the local ethics committee. We performed hotspot detection and a cortical mapping for each of four stimulation conditions: biphasic latero-medial (BI-90), biphasic posterolateral-anteromedial (BI-45), monophasic latero-medial (MO-90), and monophasic posterolateral-anteromedial (MO-45). Due to the duration of each mapping, recordings were performed in two sessions. Per session, two conditions were assessed. The order of conditions was randomized within each subject.

_Transcranial Magnetic Stimulation_

In each session, TMS was applied with a MagPro X100 with MagOption and a MCF-B70 figure-of-eight coil of 97mm diameter (Magventure, Denmark). Stimulation was performed based on a template headmodel using the (TMS Navigator, Localite, Germany). We recorded MEPs at the extensor digitorum of the right hand with surface electromyography. The muscle was located by palpation during extension of the wrist and anatomical landmarks. The skin was cleaned using 95% ethanol and abrasive gel. We used self-adhesive electrodes (Neuroline, Ambu, Germany). During the whole experiment, the subject was seated in a comfortable reclining chair and was told to relax his or her muscles.

_Hotspot Detection_

The individual hotspot for EDC representation in the right primary motor cortex was determined for each of the four conditions prior to every mapping. Based on the individual´s head anatomy, we started stimulation at the approximate location of the hotspot. Initial stimulation intensity was set to 40% of maximum stimulator output (MSO) for biphasic stimulation and 65% of MSO for monophasic stimulation [2]. If no MEPs could be elicited, intensity was increased in 5% steps. Coil position was adapted in a random fashion, while coil orientation was fixed based on the respective condition. Finally, the hotspot of each condition was defined as the spot eliciting the highest MEP with shortest latency for a given stimulation intensity.

_Resting Motor Threshold_

At this spot the resting motor threshold (RMT) was determined with the relative frequency method [3], i.e. RMT was defined as the stimulator intensity at which 5 out of 10 stimuli would elicit an MEP with an amplitude larger than 50 µVpp.

_Mapping Grid_

For mapping, a grid was created via Localite with its center 1cm anterior of the previously detected hotspot. In steps of 0.5 cm, we set 7 x 15 grid points. This resulted in a 3 cm wide grid spanning 4.5 cm anterior and 2.5 cm posterior to the hotspot. At each of the 105 grid points, 3 stimuli were ap-plied. This resulted in an average of 15 Stimuli/cm2\. The mapping was carried out for each stimula-tion condition with 110% of RMT.

_Input-Output-Curves_

After mapping, we selected the point in the grid with the highest MEP averaged over 3 stimuli (M1 spot) and identified the most anterior point eliciting any MEP (NPMA spot). To assess differences in the respective neuronal structures, we measured the stimulus-response-curve at each spot. We stimulated with 7 intensities ranging from 90% - 150% of the respective RMT in steps of 10%. We delivered 10 stimuli per intensity. The order of stimulation intensities was randomized within each subject.

_Signal Processing_

For all MEPs, latency and amplitude were estimated automatically offline, using costum-written MatLab functions. Additionally, all trials were visually inspected for correctness of the estimation. Incorrect estimates were corrected, and artifacted trials were rejected. Latency and amplitude were averaged for each grid point and for each intensity level. Additionally, we estimated for each grid point and intensity level the MEP probability based on an amplitude threshold at 50µVpp. For input-output-curve assessment, we also processed the time-course of MEPs from 5 to 60ms after the TMS pulse. The raw time-course was detrended and baselined in reference to the period from 5 to 17 ms after the TMS pulse.

_Statistical Analysis_

Regarding the input-output-curve, we calculated the influence of the categorical factors coil ORIEN-TATION (90° vs. 45°) and stimulus WAVEFORM (biphasic vs. monophasic) as well as TARGET (M1 vs. NPMA spot) on latency and amplitude, accounting for interactions between the three factors and SUBJECT as a random factor. This test was performed for each stimulus-intensity. For the time-course, we performed this analysis additionally for every time-point. Regarding the mapping, we calculated the influence of the categorical factors coil ORIENTATION and stimulus WAVEFORM on latency and amplitude, accounting for interactions between the two factors and SUBJECT as a ran-dom factor. This test was performed for each grid-point. The statistical significance of the influence was estimated by contrasting the model coefficients with a permutation analysis using 1000 repeti-tions. For mapping, we additionally performed a cluster-based permutation analysis based on the sum of coefficients of neighboring significant grid points. The significance threshold was set to 5% for all statistical tests.
