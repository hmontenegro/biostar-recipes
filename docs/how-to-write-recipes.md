# Biostar Engine


## Where to start ?

The first recommended step is to create a project where you will be given full access. To do this, click the ```Create Project``` button found at the top of your projects list. 

![alt text](https://github.com/Natay/biostar-recipes/blob/master/docs/images/create_project.png "Create Project")

This will open a page allowing you to name your project, give it an image, etc. 

Clicking `Create` makes an empty project that you have `OWNER ACCESS` to. 

Now that you have complete access to a project you can populate it with data, collaborators, and recipes. 


## 1. Add collaborators

To add collaborators , click ```Manage Access``` found on the bottom of your data list page.

![alt text](https://github.com/Natay/biostar-recipes/blob/master/docs/images/manage_access.png "Manage Access")

That will take you to an access page that allows you to add collaborators by searching for them.

![alt text](https://github.com/Natay/biostar-recipes/blob/master/docs/images/access_interface.png "Manage Access")

Only the project creator can grant `READ ACCESS` or `WRITE ACCESS` to other users.

* `READ ACCESS` Allows users to view and copy what they desire.  

* `WRITE ACCESS` Allows users to: import data, create and run recipes, and delete/edit what they create.

* `OWNER ACCESS` Given to user when creating a project. Allows them to delete/edit anything and revoke/grant access 


## 2. Import data


* Uploading files
   
   You can upload data by clicking the ```Upload File``` button found at the top of your data list page. 
   
   ![alt text](https://github.com/Natay/biostar-recipes/blob/master/docs/images/data_dash.png "Manage Access")
   
   This will open a page allowing you to name your data, give it a type, etc. 

   There is a 25 MB file size limit and every user gets approximately 300 MB of total upload space. 
   
* Copy and Paste existing data ( **Recommended** )

   Copy files or data from other projects you have `READ ACCESS` to and paste them with no limits on size.
   
   You can copy data by following these steps:
   
     1. Click on `Browse Files` to see the copying interface
     
     2. Select the files and click `Copy Data`
     
     ![alt text](https://github.com/Natay/biostar-recipes/blob/master/docs/images/copy_data.png "Copy Data")

   You can edit the information like name, type, and description by clicking the `Edit Info` button found on the bottom of every data page.
   
    ![alt text](https://github.com/Natay/biostar-recipes/blob/master/docs/images/data_info.png "Edit Data")
   
   
**Data Types**


Data types are used to sub-select for data during analysis. You can specifiy a data type when uploading or afterwards by clicking `Edit Info` . 

To change a type, you simply have to enter a string in the box labeled `Data Type`. Leaving it blank on upload will give the data a default type `DATA`.  

![alt text](https://github.com/Natay/biostar-recipes/blob/master/docs/images/data_type.png "Edit Data")

Data can have multiple types, all comma seperated. This denotes that data will be sub-selected for any type present in the comma seperated string.   

![alt text](https://github.com/Natay/biostar-recipes/blob/master/docs/images/data_type2.png "Edit Data")


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




