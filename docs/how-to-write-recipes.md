# Biostar Engine


## 3. Make your own recipes

If you have `WRITE ACCESS` or higher then you have the ability to create a recipe in one of two ways.

* Start from scratch

    To create a brand new recipe click on the ```Create Recipe``` button found at the bottom
    of your recipes list. 
    
    ![alt text](https://github.com/Natay/biostar-recipes/blob/master/docs/images/recipe_create.png "Create Recipe")
    
    This opens a page that allowing you to name your recipe, give it a picture, etc. 
    
    Clicking `Save` will create a recipe with an empty template and json file, which you can edit by clicking `View Code`. 
    
   
* Copy and Paste existing recipes ( **Recommended** )

    You can copy data by following these steps:
    
     1. Go to a recipe of choice and right under the run button should be a `Copy Recipe` button.
     
     ![alt text](https://github.com/Natay/biostar-recipes/blob/master/docs/images/recipe_copy.png "Copy Recipe")
     
     2. Click `Copy Recipe` and save it into your clipboard. 
      
     ![alt text](https://github.com/Natay/biostar-recipes/blob/master/docs/images/recipe_paste.png "Paste Recipe")
      
    This is the simplest way of populating your project with recipes.
    
    You only need `READ ACCESS` to be able to copy a recipe, so public project like `The Biostar Cookbook` are treasure troves.
    
    After pasting the recipe you can edit the code by clicking `View Code`.
    
   
**Editing recipe code**

Editing recipe code is done by clicking `View Code`. 

   * Editing json and **sub-selecting for data types**
      
      Biostar-Engine knows to look for data if `source : PROJECT`. Furthermore, `display: DROPDOWN` needs to be set for the dropdown interface correctly show.
      
      Here are simple examples that show how to sub-select:
      
      
          # Only show FASTA type in the dropdown
          { 
            name : Test
            source : PROJECT
            type : FASTA
            display: DROPDOWN
          }
          
          # Shows all data in a project
          { 
            name : Test
            source : PROJECT
            display: DROPDOWN
          }
    
         
   * Editing template
   
       Editing a template will result in the changes having to be reviewed and authorized by a staff member.
       
       For security reasons,a recipe with a changed template can not be executed without a staff member authorizing it. 
       
   ![alt text](https://github.com/Natay/biostar-recipes/blob/master/docs/images/recipe_code.png "Recipe code")
   ![alt text](https://github.com/Natay/biostar-recipes/blob/master/docs/images/recipe_code2.png "Recipe code")



## 4. Running recipes

After editing your recipe, you can finally use it to analyze some data. To do this, choose a recipe and click on the green ```Run``` button at the top of the page. 

This opens the interface that allows you to specify parameters and execute a recipe.

![alt text](https://github.com/Natay/biostar-recipes/blob/master/docs/images/run_interface.png "interface")

Clicking `Run` on an interface page starts a job in a `Queued` state

![alt text](https://github.com/Natay/biostar-recipes/blob/master/docs/images/params.png)


* Authorization requirements

     1. You need `WRITE ACCESS` to the project
     2. Recipe needs to be reviewd by a staff if any changes have been made.
     
    ![alt text](https://github.com/Natay/biostar-recipes/blob/master/docs/images/template_change.png "Recipe code")
     

* Job states

     `Queued`   | Job is staged to run.
     
     `Spooled`   | Spooling directory has been made and script is ready to be ran.
   
     `Running`   | Job script is ran using bash.
   
     `Completed` | Job is successfully completed.
   
     `Error`     | Job is unsuccessfully completed.
   
     `Deleted`   | Deleted job not appearing in your list but your recycle bin.
   
     `Paused`    | Job is paused from its previous state. 
   
     `Restored`  | A once deleted job restored from the recycle bin. 

* Deleting and Editing jobs
   
   Deleting and editing jobs are actions only allowed to the creator of the job or project. 
   
   ![alt text](https://github.com/Natay/biostar-recipes/blob/master/docs/images/job_edit.png "Recipe code")

* Gathering results

   In the same way you copy data, you can also copy result files 


## 5. Simple Demo


1. Create a project


2. Copy some data


3. Copy a recipe


4. Run recipe


5. Copy resulting files




