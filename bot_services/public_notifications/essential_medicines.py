#!/usr/bin/env python3
#
# GlobalChem - World Health Organization Bot
#
# ------------------------------------------

# Imports
# -------

import os
import json
import requests
import pandas as pd

from bs4 import BeautifulSoup

class EssentialMedicinesBot(object):

    __reference__ = [
        'Halothane', 'Isoflurane', 'Nitrous oxide (medication)', 'Oxygen therapy', 'Ketamine', 'Propofol',
        'Bupivacaine', 'Lidocaine', 'Lidocaine/epinephrine', 'Lidocaine', 'Epinephrine (medication)', 'Ephedrine',
        'Atropine', 'Midazolam', 'Morphine', 'Oxygen therapy', 'Enlarge', 'Skeletal model', 'Acetylsalicylic acid',
        'Ibuprofen', 'Paracetamol', 'Codeine', 'Fentanyl', 'Morphine', 'Methadone', 'Amitriptyline', 'Cyclizine',
        'Dexamethasone', 'Diazepam', 'Docusate sodium', 'Fluoxetine', 'Haloperidol', 'Hyoscine butylbromide',
        'Hyoscine hydrobromide', 'Lactulose', 'Loperamide', 'Metoclopramide', 'Midazolam', 'Ondansetron',
        'Senna glycosides', 'Dexamethasone', 'Epinephrine (medication)', 'Hydrocortisone', 'Loratadine',
        'Prednisolone', 'Activated charcoal (medication)', 'Acetylcysteine', 'Atropine', 'Calcium gluconate',
        'Methylene blue', 'Naloxone', 'Penicillamine', 'Prussian blue (medical use)', 'Sodium nitrite (medical use)',
        'Sodium thiosulfate (medical use)', 'Deferoxamine', 'Dimercaprol', 'Fomepizole', 'Sodium calcium edetate',
        'Succimer', 'Carbamazepine', 'Diazepam', 'Lamotrigine', 'Lorazepam', 'Magnesium sulfate (medical use)',
        'Midazolam', 'Phenobarbital', 'Phenytoin', 'Valproic acid', 'Ethosuximide', 'Valproic acid', 'Enlarge',
        'Albendazole', 'Ivermectin', 'Levamisole', 'Mebendazole', 'Niclosamide', 'Praziquantel', 'Pyrantel',
        'Albendazole', 'Diethylcarbamazine', 'Ivermectin', 'Praziquantel', 'Triclabendazole', 'Oxamniquine',
        'Albendazole', 'Mebendazole', 'Praziquantel', 'Amikacin', 'Amoxicillin', 'Amoxicillin/clavulanic acid',
        'Amoxicillin', 'Clavulanic acid', 'Ampicillin', 'Benzathine benzylpenicillin', 'Benzylpenicillin',
        'Cefalexin', 'Cefazolin', 'Chloramphenicol', 'Clindamycin', 'Cloxacillin', 'Doxycycline', 'Gentamicin',
        'Metronidazole', 'Nitrofurantoin', 'Phenoxymethylpenicillin', 'Procaine benzylpenicillin', 'Spectinomycin',
        'Sulfamethoxazole/trimethoprim', 'Sulfamethoxazole', 'Trimethoprim', 'Trimethoprim', 'Azithromycin', 'Cefixime',
        'Cefotaxime', 'Ceftriaxone', 'Cefuroxime', 'Ciprofloxacin', 'Clarithromycin', 'Piperacillin/tazobactam',
        'Piperacillin', 'Tazobactam', 'Vancomycin', 'Ceftazidime', 'Meropenem', 'Vancomycin', 'Cefiderocol',
        'Ceftazidime/avibactam', 'Ceftazidime', 'Avibactam', 'Colistin', 'Fosfomycin', 'Linezolid', 'Meropenem/vaborbactam',
        'Meropenem', 'Vaborbactam', 'Plazomicin', 'Polymyxin B', 'Clofazimine', 'Dapsone', 'Rifampicin', 'Enlarge',
        'Ethambutol', 'Ethambutol/isoniazid/pyrazinamide/rifampicin', 'Ethambutol', 'Isoniazid', 'Pyrazinamide',
        'Rifampicin', 'Ethambutol/isoniazid/rifampicin', 'Ethambutol', 'Isoniazid', 'Rifampicin', 'Isoniazid',
        'Isoniazid/pyrazinamide/rifampicin', 'Isoniazid', 'Pyrazinamide', 'Rifampicin', 'Isoniazid/rifampicin',
        'Isoniazid', 'Rifampicin', 'Isoniazid/rifapentine (page does not exist)', 'Isoniazid', 'Rifapentine',
        'Moxifloxacin', 'Pyrazinamide', 'Rifabutin', 'Rifampicin', 'Rifapentine', 'Amikacin', 'Amoxicillin/clavulanic acid',
        'Amoxicillin', 'Clavulanic acid', 'Bedaquiline', 'Clofazimine', 'Cycloserine', 'Delamanid', 'Ethionamide',
        'Levofloxacin', 'Linezolid', 'Meropenem', 'Moxifloxacin', 'P-aminosalicylic acid', 'Streptomycin', 'Amphotericin B',
        'Clotrimazole', 'Fluconazole', 'Flucytosine', 'Griseofulvin', 'Itraconazole', 'Nystatin', 'Voriconazole', 'Micafungin',
        'Potassium iodide', 'Aciclovir', 'Abacavir', 'Lamivudine', 'Tenofovir disoproxil fumarate', 'Zidovudine', 'Efavirenz',
        'Nevirapine', 'Enlarge', 'Atazanavir/ritonavir', 'Atazanavir', 'Ritonavir', 'Darunavir', 'Lopinavir/ritonavir',
        'Lopinavir', 'Ritonavir', 'Ritonavir', 'Dolutegravir', 'Raltegravir', 'Abacavir/lamivudine', 'Abacavir', 'Lamivudine',
        'Dolutegravir/lamivudine/tenofovir', 'Dolutegravir', 'Lamivudine', 'Tenofovir', 'Efavirenz/emtricitabine/tenofovir',
        'Efavirenz/lamivudine/tenofovir', 'Efavirenz', 'Lamivudine', 'Tenofovir', 'Emtricitabine/tenofovir', 'Emtricitabine',
        'Tenofovir', 'Lamivudine/zidovudine', 'Lamivudine', 'Zidovudine', 'Isoniazid/pyridoxine/sulfamethoxazole/trimethoprim',
        'Isoniazid', 'Pyridoxine', 'Sulfamethoxazole', 'Trimethoprim', 'Ribavirin', 'Valganciclovir', 'Oseltamivir',
        'Valganciclovir', 'Entecavir', 'Tenofovir disoproxil fumarate', 'Daclatasvir', 'Daclatasvir/sofosbuvir',
        'Daclatasvir', 'Sofosbuvir', 'Glecaprevir/pibrentasvir', 'Glecaprevir', 'Pibrentasvir', 'Sofosbuvir',
        'Sofosbuvir/velpatasvir', 'Sofosbuvir', 'Velpatasvir', 'Dasabuvir', 'Ledipasvir/sofosbuvir', 'Ledipasvir',
        'Sofosbuvir', 'Ombitasvir/paritaprevir/ritonavir', 'Ombitasvir', 'Paritaprevir', 'Ritonavir', 'Ribavirin',
        'Pegylated interferon-alpha-2a', 'Pegylated interferon-alpha-2b', 'Diloxanide', 'Metronidazole', 'Amphotericin B',
        'Miltefosine', 'Paromomycin', 'Sodium stibogluconate', 'Meglumine antimoniate', 'Amodiaquine', 'Artemether',
        'Artemether/lumefantrine', 'Artemether', 'Lumefantrine', 'Artesunate', 'Artesunate/amodiaquine', 'Artesunate',
        'Amodiaquine', 'Artesunate/mefloquine', 'Artesunate', 'Mefloquine', 'Artesunate/pyronaridine tetraphosphate',
        'Artesunate', 'Pyronaridine tetraphosphate', 'Chloroquine', 'Dihydroartemisinin/piperaquine phosphate',
        'Dihydroartemisinin', 'Piperaquine phosphate', 'Doxycycline', 'Mefloquine', 'Primaquine', 'Quinine',
        'Sulfadoxine/pyrimethamine', 'Sulfadoxine', 'Pyrimethamine', 'Amodiaquine', 'Sulfadoxine/pyrimethamine',
        'Chloroquine', 'Doxycycline', 'Mefloquine', 'Proguanil', 'Sulfadoxine/pyrimethamine', 'Sulfadoxine',
        'Pyrimethamine', 'Pyrimethamine', 'Sulfadiazine', 'Sulfamethoxazole/trimethoprim', 'Sulfamethoxazole',
        'Trimethoprim', 'Pentamidine', 'Fexinidazole', 'Pentamidine', 'Suramin sodium', 'Eflornithine', 'Melarsoprol',
        'Nifurtimox', 'Melarsoprol', 'Benznidazole', 'Nifurtimox', 'Ivermectin', 'Acetylsalicylic acid', 'Ibuprofen',
        'Paracetamol', 'Sumatriptan', 'Propranolol', 'Adalimumab', 'Azathioprine', 'Ciclosporin', 'Tacrolimus',
        'Arsenic trioxide', 'Asparaginase', 'Bendamustine', 'Bleomycin', 'Calcium folinate', 'Capecitabine', 'Carboplatin',
        'Chlorambucil', 'Cisplatin', 'Cyclophosphamide', 'Cytarabine', 'Dacarbazine', 'Dactinomycin', 'Daunorubicin',
        'Docetaxel', 'Doxorubicin', 'Etoposide', 'Fludarabine', 'Fluorouracil', 'Gemcitabine', 'Hydroxycarbamide',
        'Ifosfamide', 'Irinotecan', 'Melphalan', 'Mercaptopurine', 'Methotrexate', 'Oxaliplatin', 'Paclitaxel', 'Pegaspargase',
        'Procarbazine', 'Realgar-Indigo naturalis', 'Tioguanine', 'Vinblastine', 'Vincristine', 'Vinorelbine', 'Tretinoin',
        'Tretinoin', 'Bortezomib', 'Dasatinib', 'Erlotinib', 'Everolimus', 'Ibrutinib', 'Imatinib', 'Nilotinib', 'Rituximab',
        'Trastuzumab', 'Filgrastim', 'Lenalidomide', 'Nivolumab', 'Thalidomide', 'Abiraterone', 'Anastrozole', 'Bicalutamide',
        'Dexamethasone', 'Hydrocortisone', 'Leuprorelin', 'Methylprednisolone', 'Prednisolone', 'Tamoxifen', 'Allopurinol',
        'Mesna', 'Rasburicase', 'Zoledronic acid', 'Biperiden', 'Levodopa/carbidopa', 'Levodopa', 'Carbidopa', 'Iron supplement',
        'Ferrous salt/folic acid', 'Ferrous salt', 'Folic acid', 'Folic acid', 'Hydroxocobalamin', 'Erythropoiesis-stimulating agent',
        'Dabigatran', 'Enoxaparin', 'Heparin sodium', 'Phytomenadione', 'Protamine sulfate', 'Tranexamic acid', 'Warfarin',
        'Desmopressin', 'Heparin sodium', 'Protamine sulfate', 'Warfarin', 'Deferoxamine', 'Hydroxycarbamide', 'Enlarge',
        'Fresh frozen plasma', 'Platelet concentrates', 'Packed red blood cell', 'Whole blood', 'Rho(D) immune globulin',
        'Anti-rabies immunoglobulin', 'Anti-tetanus immunoglobulin', 'Human normal immunoglobulin', 'Factor VIII (medication)',
        'Factor IX complex', 'Dextran 70', 'Bisoprolol', 'Glyceryl trinitrate (pharmacology)', 'Isosorbide dinitrate',
        'Verapamil', 'Bisoprolol', 'Digoxin', 'Epinephrine (medication)', 'Lidocaine', 'Verapamil', 'Amiodarone', 'Amlodipine',
        'Bisoprolol', 'Enalapril', 'Hydralazine', 'Hydrochlorothiazide', 'Lisinopril/amlodipine', 'Lisinopril', 'Amlodipine',
        'Lisinopril/hydrochlorothiazide', 'Lisinopril', 'Hydrochlorothiazide', 'Losartan', 'Methyldopa', 'Telmisartan/amlodipine',
        'Telmisartan', 'Amlodipine', 'Telmisartan/hydrochlorothiazide', 'Telmisartan', 'Hydrochlorothiazide', 'Sodium nitroprusside',
        'Bisoprolol', 'Digoxin', 'Enalapril', 'Furosemide', 'Hydrochlorothiazide', 'Losartan', 'Spironolactone',
        'Dopamine (medication)', 'Acetylsalicylic acid', 'Clopidogrel', 'Alteplase', 'Streptokinase', 'Simvastatin',
        'Miconazole', 'Selenium disulfide', 'Sodium thiosulfate (medical use)', 'Terbinafine', 'Mupirocin',
        'Potassium permanganate (medical use)', 'Silver sulfadiazine', 'Betamethasone', 'Calamine', 'Hydrocortisone',
        'Benzoyl peroxide', 'Calcipotriol', 'Coal tar', 'Fluorouracil', 'Podophyllum resin', 'Salicylic acid (medical use)',
        'Urea-containing cream', 'Benzyl benzoate', 'Permethrin', 'Fluorescein (medical use)', 'Tropicamide', 'Amidotrizoate',
        'Barium sulfate suspension', 'Iohexol', 'Barium sulfate suspension', 'Meglumine iotroxate', 'Chlorhexidine',
        'Alcohol (medical use)', 'Povidone-iodine', 'Alcohol based hand rub', 'Chlorine base compound', 'Chloroxylenol',
        'Glutaral', 'Amiloride', 'Furosemide', 'Hydrochlorothiazide', 'Mannitol', 'Spironolactone', 'Hydrochlorothiazide',
        'Mannitol', 'Spironolactone', 'Pancreatic enzymes (medication)', 'Omeprazole', 'Ranitidine', 'Dexamethasone',
        'Metoclopramide', 'Ondansetron', 'Aprepitant', 'Sulfasalazine', 'Hydrocortisone', 'Prednisolone', 'Senna glycosides',
        'Oral rehydration salts', 'Zinc sulfate (medical use)', 'Oral rehydration salts', 'Zinc sulfate (medical use)',
        'Fludrocortisone', 'Hydrocortisone', 'Testosterone (medication)', 'Medroxyprogesterone acetate', 'Insulin injection (soluble)',
        'Intermediate-acting insulin', 'Long-acting insulin analogues', 'Empagliflozin', 'Gliclazide', 'Metformin', 'Metformin',
        'Glucagon (medication)', 'Diazoxide', 'Levothyroxine', 'Potassium iodide', 'Methimazole', 'Propylthiouracil', "Lugol's solution", 'Methimazole', 'Potassium iodide',
        'Propylthiouracil', 'Tuberculin', 'Antivenom immunoglobulin', 'Diphtheria antitoxin', 'Rabies immunoglobulin', 'Enlarge', 'Bacillus Calmette–Guérin', 'Diphtheria vaccine', 'Haemophilus influenzae type b vaccine', 'Hepatitis B vaccine',
        'HPV vaccine', 'Measles vaccine', 'Pertussis vaccine', 'Pneumococcal vaccine', 'Poliomyelitis vaccine', 'Rotavirus vaccine', 'Rubella vaccine', 'Tetanus vaccine', 'Japanese encephalitis vaccine', 'Tick-borne encephalitis vaccine', 'Yellow fever vaccine', 'Cholera vaccine', 'Dengue vaccine', 'Hepatitis A vaccine',
        'Meningococcal meningitis vaccine', 'Rabies vaccine', 'Typhoid vaccine', 'Influenza vaccine', 'Mumps vaccine', 'Varicella vaccine', 'Atracurium', 'Neostigmine', 'Suxamethonium', 'Vecuronium', 'Pyridostigmine', 'Vecuronium', 'Aciclovir', 'Azithromycin', 'Erythromycin', 'Gentamicin', 'Natamycin',
        'Ofloxacin', 'Tetracycline', 'Prednisolone', 'Tetracaine', 'Acetazolamide', 'Latanoprost', 'Pilocarpine', 'Timolol', 'Atropine', 'Epinephrine (medication)', 'Bevacizumab', 'Ethinylestradiol/levonorgestrel', 'Ethinylestradiol', 'Levonorgestrel', 'Ethinylestradiol/norethisterone', 'Ethinylestradiol', 'Norethisterone', 'Levonorgestrel', 'Ulipristal', 'Estradiol cypionate/medroxyprogesterone acetate', 'Estradiol cypionate', 'Medroxyprogesterone acetate', 'Medroxyprogesterone acetate', 'Norethisterone enantate', 'IUD with copper', 'IUD with progestogen', 'Condom', 'Diaphragm (contraceptive)', 'Etonogestrel birth control implant', 'Levonorgestrel-releasing implant', 'Ethinylestradiol/etonogestrel', 'Ethinylestradiol', 'Etonogestrel', 'Progesterone vaginal ring', 'Clomifene', 'Carbetocin', 'Ergometrine', 'Mifepristone', 'Misoprostol', 'Misoprostol', 'Oxytocin (medication)', 'Nifedipine', 'Dexamethasone', 'Multivitamin', 'Tranexamic acid', 'Caffeine citrate', 'Chlorhexidine', 'Ibuprofen', 'Prostaglandin E1', 'Pulmonary surfactant (medication)', 'Intraperitoneal dialysis solution', 'Chlorpromazine', 'Fluphenazine', 'Haloperidol', 'Paliperidone', 'Risperidone', 'Chlorpromazine', 'Clozapine', 'Haloperidol', 'Amitriptyline', 'Fluoxetine', 'Fluoxetine', 'Sertraline', 'Carbamazepine', 'Lithium (medication)', 'Valproic acid', 'Diazepam', 'Clomipramine', 'Bupropion', 'Nicotine replacement therapy', 'Varenicline', 'Methadone', 'Budesonide', 'Budesonide/formoterol', 'Budesonide', 'Formoterol', 'Epinephrine (medication)', 'Ipratropium bromide', 'Salbutamol', 'Tiotropium', 'Oral rehydration salts', 'Potassium chloride (medical use)', 'Intravenous sugar solution', 'Intravenous sugar solution', 'Potassium chloride (medical use)', 'Saline (medicine)', 'Intravenous sodium bicarbonate', 'Sodium lactate, compound solution', 'Water for injection', 'Ascorbic acid', 'Calcium supplement', 'Colecalciferol', 'Ergocalciferol', 'Iodine (medical use)', 'Multiple micronutrient powder', 'Nicotinamide', 'Pyridoxine', 'Retinol', 'Riboflavin', 'Thiamine', 'Calcium gluconate', 'Acetic acid (medical use)', 'Budesonide', 'Ciprofloxacin', 'Xylometazoline', 'Allopurinol', 'Chloroquine', 'Azathioprine', 'Hydroxychloroquine', 'Methotrexate', 'Penicillamine', 'Sulfasalazine', 'Aspirin', 'Fluoride therapy', 'Glass ionomer cement', 'Silver diamine fluoride', 'Cost-benefit ratio', 'Thiopental', 'Hydromorphone', 'Oxycodone', 'Dolasetron', 'Granisetron', 'Palonosetron', 'Tropisetron', 'Cetirizine', 'Fexofenadine', 'Prednisone', 'Diazepam', 'Midazolam', 'Erythromycin', 'Imipenem/cilastatin', 'Terizidone', 'Prothionamide', 'Aspergillosis', 'Histoplasmosis', 'Sporotrichosis', 'Paracoccidioidomycosis', 'Talaromyces marneffei', 'Chromoblastomycosis', 'Anidulafungin', 'Caspofungin', 'Valaciclovir', 'Cytomegalovirus', 'Sofosbuvir', 'Daclatasvir', 'Ribavirin', 'Tinidazole', 'Artesunate', 'Amodiaquine', 'Mefloquine', 'Sulfadoxine', 'Pyrimethamine', 'Certolizumab pegol', 'Etanercept', 'Golimumab', 'Infliximab', 'Biosimilar', 'Afatinib', 'Gefitinib', 'Pembrolizumab', 'Biosimilar', 'Enzalutamide', 'Flutamide', 'Nilutamide', 'Goserelin', 'Triptorelin', 'Prednisone', 'Trihexyphenidyl', 'Benserazide', 'Epoetin alfa', 'Epoetin beta', 'Epoetin theta', 'Darbepoetin alfa', 'Methoxy polyethylene glycol-epoetin beta', 'Apixaban', 'Edoxaban', 'Rivaroxaban', 'Dalteparin', 'Nadroparin', 'Deferasirox', 'Carvedilol', 'Metoprolol', 'Atenolol', 'Carvedilol', 'Metoprolol', 'Chlorothiazide', 'Chlorthalidone', 'Indapamide', 'Bumetanide', 'Torasemide', 'Atorvastatin', 'Fluvastatin', 'Lovastatin', 'Pravastatin', 'Calcitriol', 'Tacalcitol', 'Podophyllotoxin', 'Atropine', 'Cyclopentolate', 'Propanol', 'Iodine', 'Bumetanide', 'Torasemide', 'Chlorothiazide', 'Chlorthalidone', 'Mesalazine', 'Bisacodyl', 'Norethisterone', 'Insulin degludec', 'Insulin detemir', 'Insulin glargine', 'Canagliflozin', 'Dapagliflozin', 'Carbimazole', 'Amikacin', 'Kanamycin', 'Netilmicin', 'Tobramycin', 'Chlortetracycline', 'Oxytetracycline', 'Carbachol', 'Cyclopentolate hydrochloride', 'Homatropine hydrobromide', 'Methylergometrine', 'Abortion law', 'Abortion debate', 'Indometacin', 'Prostaglandin E2', 'Risperidone', 'Citalopram', 'Escitalopram', 'Fluvoxamine', 'Paroxetine', 'Sertraline', 'Buprenorphine', 'Beclometasone', 'Ciclesonide', 'Flunisolide', 'Fluticasone', 'Mometasone', 'Beclometasone/formoterol', 'Budesonide/salmeterol (page does not exist)', 'Fluticasone/formoterol (page does not exist)',
        'Fluticasone furoate/vilanterol', 'Mometasone/formoterol', 'Terbutaline', 'Aclidinium', 'Glycopyrronium', 'Umeclidinium', 'Ergocalciferol', 'Colecalciferol', 'Ofloxacin'
    ]

    def __init__(self):

        self.name = 'essential_medicines'
        self.master_url = "https://en.wikipedia.org/wiki/WHO_Model_List_of_Essential_Medicines"

        self.essential_medicine_lists = []
        self.global_chem_lists = []

    def get_essential_medicines_list(self):

        '''

        Get the Essential Medicines List

        '''

        page = requests.get(self.master_url)
        soup = BeautifulSoup(page.content, "html.parser")

        essential_medicines = soup.find_all("a")

        start = False
        molecules = []

        for essential_medicine in essential_medicines:

            title = essential_medicine.get('title', 'N/A')
            if title == 'N/A':
                continue
            elif 'edit section: references' in title.lower():
                break
            elif 'edit section: inhalational medicines' in title.lower():
                start = True
                continue
            elif not start:
                continue
            elif 'edit section' in title.lower():
                continue
            elif len(title) < 3:
                continue


            molecules.append(title)

        self.essential_medicine_lists = molecules


    def check_list_status(self):

        '''

        File a Github Issue on the global chem Repo if there is a misalignment.

        '''

        differences = set(self.essential_medicine_lists).symmetric_difference(set(self.__reference__))

        if len(differences) > 1:

            token = os.getenv('GITHUB_TOKEN')

            headers = {
                "Accept": 'Accept: application/vnd.github+json',
                "Authorization" : "token {}".format(token)
            }

            data = {
                "title": "Public Notification System: Essential Medicines from WHO",
                "body": "Differences found that need to be inspected: %s" % differences

            }

            username = 'Sulstice'
            repository_name = 'global-chem'

            url = "https://api.github.com/repos/{}/{}/issues".format(
                username,
                repository_name
            )

            response = requests.post(
                url,
                data=json.dumps(data),
                headers=headers
            )


if __name__ == '__main__':

    bot = EssentialMedicinesBot()
    bot.get_essential_medicines_list()
    bot.check_list_status()