### **HackBio Internship: Stage 3**

**Authors (@slack):**   
Lakshana Bakthavachalam **(@Lakshana),**  Nada ElSayed Ahmed **(@Nada\_EA)**

**GitHub link to markdown  file: [https://github.com/Nada-EA/HackBio2024-AMR/blob/main/Stage3/Final%20Report.md]**

**GitHub link to Stage 3 repository: [https://github.com/Nada-EA/HackBio2024-AMR/tree/main/Stage3]**

**Phase I**

**Google Doc:** [https://docs.google.com/document/d/1V0sQZQAMyaK9svz\_VcVEsrdnY6KWJLA7qIHcCffQY0A/edit\#heading=h.6o1oyespi2s2](https://docs.google.com/document/d/1V0sQZQAMyaK9svz_VcVEsrdnY6KWJLA7qIHcCffQY0A/edit#heading=h.6o1oyespi2s2)

**Phase II**

**Global Cholera Outbreaks (1949 \- 2016\)**

The “Cholera outbreak” data was obtained from the Global Health Observatory (GHO), a data repository of WHO, which provides information about the health-related statistics for its 194 member countries. The data was pre-processed and visualized using R, based on three key parameters: number of cases, number of deaths, and case fatality rate.


 ![map_cases](https://github.com/user-attachments/assets/cc93b0e1-e6fe-4b1e-85cc-d2e5af7dfa69)
 
**Fig1: HOTSPOTS OF CHOLERA OUTBREAK**


Between 1949 and 2016, India had the largest cholera outbreaks, followed by Peru and Democratic Republic of Congo. Other countries like Afghanistan, Bangladesh,  Malawi, Madagascar, Somalia, Tanzania, also experienced an outbreak but it was not alarming (in terms of scale) when compared to India and Peru. Among the listed countries, Bangladesh and DRC are classified as endemic regions.


![Country_total_metric](https://github.com/user-attachments/assets/dfef7f4b-ef5a-4540-a781-f59379f52af8)

**Fig2: TOP 20 COUNTRIES MOST AFFECTED BY CHOLERA**


The bar graphs highlight the top 20 countries that were most affected by cholera. Despite the alarming case numbers, the fatality rate is found to be under control. For instance, the number of cases exceeds over 1 million in India; however, the fatality is less than 50% of the case total. A similar pattern is also observed in Bangladesh.


![Cases_years](https://github.com/user-attachments/assets/c54e2695-418b-49d5-ac90-14be78468f77)

**Fig3: TRENDS IN CHOLERA CASES \- 1949 TO 2016**


In 1949, approximately 0.2 million cases of cholera were reported worldwide. This number sharply decreased to 0.1 million by 1954 and remained as such until 1990, with minor fluctuations. However, in the year 1991, the cases surged and attained a highest peak of 0.6 million before gradually declining to around 0.1 million by the end of 1999\. From 2000 to 2010, the trend stabilized between 0.1 and 0.2 million cases, until a peak similar to that of 1991 occurred again in 2011, after which it diminished back to normal level by 2016\.

###  
 ![Case_Fatality_years](https://github.com/user-attachments/assets/a991b7ed-68e9-4121-93a9-2e9750bbcbe9)
 
### **Fig4: TRENDS IN CASE FATALITY  \- 1949 TO 2016**


From 1949 to 1960, the average case fatality rate was quite high. However, after 1960, it gradually declined and stabilized at under 10 (in number) throughout the following years. With reference to Fig3, the case number surged to 0.6 million in both 1991 and 2011, whereas the case-fatality rate during those years were low, reflecting advancements in medical development and public health interventions.



### **Cholera Outbreak Report Template**

### **Focus Country: \[Country name\] (Year range)**

### **Key Metrics**

### **Total Cases:** \[Insert number\]

### **Total Fatalities:** \[Insert number\]

### **Average Case-Fatality Rate:** \[Insert percentage\]

**Summary:**

Between \[start year\] and \[end year\], \[insert country\] experienced periodic cholera outbreaks, with particularly severe epidemics in \[insert relevant years\]. These outbreaks underscore the ongoing need for sustained efforts in prevention, early detection, and timely response. Despite fluctuations in cases and fatalities, there has been a notable improvement in case-fatality rates, suggesting advancements in medical care and outbreak management.

### **Cholera Outbreak Graph:** \[Insert graph visualizing cholera outbreak data, such as total cases, fatalities, and case-fatality rate trends over time.\]



### **Example: Focus Country: Nigeria (1970 \- 2016\)**

**Key Metrics**

**Total Cases:** 310217

**Total Fatalities:** 21479

**Average Case-Fatality Rate:** 6.92%

**Summary:**  
Between 1970 and 2016, Nigeria experienced periodic cholera outbreaks, with particularly severe epidemics in the early 1990s and early 2010s. These outbreaks underscore the ongoing need for sustained efforts in prevention, early detection, and timely response to mitigate future outbreaks. However, no data is available for Nigeria prior to 1970, which limits insight into earlier trends. Despite fluctuations in the number of cases and fatalities over the years, there has been a notable improvement in the case-fatality rate, suggesting advancements in medical care and outbreak management.


### **Cholera Outbreak Graph:** 

![Nigeria](https://github.com/user-attachments/assets/7b0fde2a-185b-422d-94d0-654c3bc34659)



**Shiny App: CholeraWatch**

**Link to app: [https://jpjxx1-nada-ahmed.shinyapps.io/cholerawatch/](https://jpjxx1-nada-ahmed.shinyapps.io/cholerawatch/)**

![CholeraWatch_app](https://github.com/user-attachments/assets/4e013d3b-712f-46d2-bb6b-e2dea721c6aa)


The CholeraWatch interactive app enables users to explore cholera outbreaks globally, based on WHO data from 1949 to 2016\. Users can select a country from a list and view key metrics, including total cases, total fatalities, and the case-fatality rate, which are displayed in three info boxes. Additionally, users can customize the time range to analyze trends over specific periods. For each selected country, the app also generates a plot with three graphs, visualizing the same metrics—total cases, total fatalities, and case-fatality rate—over time.

**Data Source:** 

World Health Organization Global Health Observatory data repository: Cholera  
[https://apps.who.int/gho/data/node.main.174?lang=en](https://apps.who.int/gho/data/node.main.174?lang=en)  
